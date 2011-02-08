#include "ffmpeg.hh"
#include "config.hh"
#include "util.hh"
#include <iostream>
#include <stdexcept>
#include <QtGlobal>

// Somehow ffmpeg headers give errors that these are not defined...
#define INT64_C Q_INT64_C
#define UINT64_C Q_UINT64_C

extern "C" {
#include AVCODEC_INCLUDE
#include AVFORMAT_INCLUDE
}

/// A custom allocator that uses av_malloc for aligned buffers
template <typename T> class AvMalloc: public std::allocator<T> {
public:
/*	pointer allocate(size_type count, allocator<void>::const_pointer* = 0) {
		pointer ptr = static_cast<pointer>(m(count * sizeof(T)));
		if (!ptr) throw std::bad_alloc();
		return ptr;
	}
    void deallocate(pointer ptr, size_type) { f(ptr); }*/
private:
	void* m(size_t n) { return av_malloc(n); }
	void f(void* ptr) { av_free(ptr); }
};


/*static*/ QMutex FFmpeg::s_avcodec_mutex;

FFmpeg::FFmpeg(std::string const& _filename):
  m_filename(_filename), m_quit(), m_running(), m_eof(),
  pFormatCtx(), pAudioCodecCtx(), pAudioCodec(),
  audioStream(-1), m_position()
{
	open(); // Throws on error
	m_running = true;
	start();
}

FFmpeg::~FFmpeg() {
	m_quit = true;
	audioQueue.setEof();
	audioQueue.reset();
	wait();
	// TODO: use RAII for freeing resources (to prevent memory leaks)
	QMutexLocker l(&s_avcodec_mutex); // avcodec_close is not thread-safe
	if (pAudioCodecCtx) avcodec_close(pAudioCodecCtx);
	if (pFormatCtx) av_close_input_file(pFormatCtx);
}

double FFmpeg::duration() const {
	double d = m_running ? pFormatCtx->duration / double(AV_TIME_BASE) : getNaN();
	return d >= 0.0 ? d : getInf();
}

void FFmpeg::open() {
	QMutexLocker l(&s_avcodec_mutex);
	av_register_all();
	av_log_set_level(AV_LOG_ERROR);
	if (av_open_input_file(&pFormatCtx, m_filename.c_str(), NULL, 0, NULL)) throw std::runtime_error("Cannot open input file");
	if (av_find_stream_info(pFormatCtx) < 0) throw std::runtime_error("Cannot find stream information");
	pFormatCtx->flags |= AVFMT_FLAG_GENPTS;
	audioStream = -1;
	// Take the first video/audio streams
	for (unsigned int i=0; i<pFormatCtx->nb_streams; i++) {
		AVCodecContext* cc = pFormatCtx->streams[i]->codec;
		cc->workaround_bugs = FF_BUG_AUTODETECT;
		if (audioStream == -1 && cc->codec_type==CODEC_TYPE_AUDIO) audioStream = i;
	}
	if (audioStream == -1) throw std::runtime_error("No audio stream found");
	AVCodecContext* cc = pFormatCtx->streams[audioStream]->codec;
	pAudioCodec = avcodec_find_decoder(cc->codec_id);
	audioQueue.setRateChannels(cc->sample_rate, cc->channels);
	if (!pAudioCodec) throw std::runtime_error("Cannot find audio codec");
	if (avcodec_open(cc, pAudioCodec) < 0) throw std::runtime_error("Cannot open audio codec");
	pAudioCodecCtx = cc;
}

void FFmpeg::run() {
	int errors = 0;
	while (!m_quit) {
		try {
			if (m_seekTarget == m_seekTarget) seek_internal();
			decodeNextFrame();
			m_eof = false;
			errors = 0;
		} catch (eof_error&) {
			audioQueue.setEof();
			m_eof = true;
			msleep(100);
		} catch (std::exception& e) {
			std::cerr << "FFMPEG error: " << e.what() << std::endl;
			if (++errors > 2) { std::cerr << "FFMPEG terminating due to errors" << std::endl; m_quit = true; }
		}
	}
	audioQueue.setEof();
	m_running = false;
	m_eof = true;
}

void FFmpeg::seek(double time, bool wait) {
	m_seekTarget = time;
	audioQueue.reset(); // Empty these to unblock the internals in case buffers were full
	if (wait) while (!m_quit && m_seekTarget == m_seekTarget) msleep(10);
}

void FFmpeg::seek_internal() {
	audioQueue.reset();
	int flags = 0;
	// FIXME: if (m_seekTarget < position()) flags |= AVSEEK_FLAG_BACKWARD;
	int stream = -1;
	stream = audioStream;
	int64_t target = m_seekTarget * AV_TIME_BASE;
	const AVRational time_base_q = { 1, AV_TIME_BASE };  // AV_TIME_BASE_Q is the same thing with C99 struct literal (not supported by MSVC)
	if (stream != -1) target = av_rescale_q(target, time_base_q, pFormatCtx->streams[stream]->time_base);
	av_seek_frame(pFormatCtx, stream, target, flags);
	m_seekTarget = getNaN(); // Signal that seeking is done
}

void FFmpeg::decodeNextFrame() {
	struct ReadFramePacket: public AVPacket {
		AVFormatContext* m_s;
		ReadFramePacket(AVFormatContext* s): m_s(s) {
			if (av_read_frame(s, this) < 0) throw eof_error();
		}
		~ReadFramePacket() { av_free_packet(this); }
		double time() {
			return uint64_t(dts) == uint64_t(AV_NOPTS_VALUE) ?
			  getNaN() : double(dts) * av_q2d(m_s->streams[stream_index]->time_base);
		}
	};

	struct AVFrameWrapper {
		AVFrame* m_frame;
		AVFrameWrapper(): m_frame(avcodec_alloc_frame()) {
			if (!m_frame) throw std::runtime_error("Unable to allocate AVFrame");
		}
		~AVFrameWrapper() { av_free(m_frame); }
		operator AVFrame*() { return m_frame; }
		AVFrame* operator->() { return m_frame; }
	} videoFrame;

	int frameFinished=0;
	while (!frameFinished) {
		ReadFramePacket packet(pFormatCtx);
		int packetSize = packet.size;
		quint8 *packetData = packet.data;
		if (packet.stream_index==audioStream) {
			while (packetSize) {
				if (packetSize < 0) throw std::logic_error("negative audio packet size?!");
				if (m_quit || m_seekTarget == m_seekTarget) return;
				std::vector<char, AvMalloc<char> > decodedBuffer(AVCODEC_MAX_AUDIO_FRAME_SIZE);
				int outsize = AVCODEC_MAX_AUDIO_FRAME_SIZE;
				int bytesUsed;
				{
					short* decodedPtr = reinterpret_cast<short*>(&decodedBuffer[0]);
#if LIBAVCODEC_VERSION_INT > ((52<<16)+(25<<8)+0)
					bytesUsed = avcodec_decode_audio3(pAudioCodecCtx, decodedPtr, &outsize, &packet);
#else
					bytesUsed = avcodec_decode_audio2(pAudioCodecCtx, decodedPtr, &outsize, packetData, packetSize);
#endif
				}
				// Negative means an error
				if (bytesUsed < 0) throw std::runtime_error("cannot decode audio frame");
				// Move forward within the packet
				packetSize -= bytesUsed;
				packetData += bytesUsed;
				// Update position if timecode is available
				if (packet.time() == packet.time()) m_position = packet.time();
				// Output samples
#define OUTPUT_SAMPLES(TYPE, MAX_VALUE) \
	{\
		outsize /= sizeof(TYPE); /* Convert bytes into samples */ \
		TYPE* samples = reinterpret_cast<TYPE*>(&decodedBuffer[0]);\
		audioQueue.input(samples, samples + outsize, 1.0 / MAX_VALUE);\
	}
				switch (pAudioCodecCtx->sample_fmt) {
					case SAMPLE_FMT_S16: OUTPUT_SAMPLES(short, 32767.0); break;
					case SAMPLE_FMT_S32: OUTPUT_SAMPLES(int, 8388607.0); break; // 24 bit samples padded to 32 bits
					case SAMPLE_FMT_FLT: OUTPUT_SAMPLES(float, 1.0); break;
					case SAMPLE_FMT_DBL: OUTPUT_SAMPLES(double, 1.0); break;
					default: throw std::runtime_error("Unsupported sample format");
				}
#undef OUTPUT_SAMPLES
				m_position += outsize / audioQueue.samplesPerSecond();  // New position in case the next packet doesn't have packet.time()
			}
			// Audio frames are always finished
			frameFinished = 1;
		}
	}
}

