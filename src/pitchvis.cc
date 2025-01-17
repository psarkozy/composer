
#include "notegraphwidget.hh"
#include "pitchvis.hh"
#include "pitch.hh"
#include "ffmpeg.hh"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <QPainter>
#include <QProgressDialog>
#include <QLabel>
#include <QSettings>
#include <QDebug>

PitchVis::PitchVis(QString const& filename, QWidget *parent, int visId)
	: QThread(parent), mutex(), fileName(filename), duration(), moreAvailable(), quit(),
	  cancelled(), restart(), m_x1(), m_y1(), m_x2(), m_y2(), m_visId(visId), condition()
{
    FFTImage = QImage();
	start(); // Launch the thread
}

void PitchVis::stop()
{
	QMutexLocker locker(&mutex);
	quit = true;
	condition.wakeOne();
}

void PitchVis::cancel()
{
	QMutexLocker locker(&mutex);
	cancelled = true;
}

void PitchVis::run()
{
	bool analyzingSuccess = false;
	try {
		// Initialize FFmpeg decoding
		std::string file(fileName.toLocal8Bit().data(), fileName.toLocal8Bit().size());
		FFmpeg mpeg(file);
		{
			QMutexLocker locker(&mutex);
			paths.clear();
			position = 0.0;
			duration = mpeg.duration(); // Estimation
		}
		unsigned rate = mpeg.audioQueue.getRate();

        //rate = 44100 ; //BECAUSE FUCK YOU THATS WHY
		unsigned channels = mpeg.audioQueue.getChannels();
		if (channels == 0) throw std::runtime_error("No audio channels found");
        qDebug() << "Loaded file has a sample rate of "<< rate <<" and a channel count of " <<channels<< " duration of "<<duration ;
		std::vector<Analyzer> analyzers(channels, Analyzer(rate, ""));
		// Process the entire song
		std::vector<float> data;
		data.reserve((duration + 1.0) * rate * channels);
		unsigned x = 0;
		while (mpeg.audioQueue.output(data)) {
			// Process as much as can be processed at this point
			while (data.size() / channels - x >= analyzers[0].processSize()) {
				// Pitch detection
				for (unsigned ch = 0; ch < channels; ++ch) {
					analyzers[ch].process(da::step_iterator<float>(&data[x * channels + ch], channels));
				}
				x += analyzers[0].processStep();
				// Update progress and check for quit flag
				QMutexLocker locker(&mutex);
				if (quit) return;
				else if (cancelled) break;
				double t = analyzers[0].getTime();
				position = t;
				duration = std::max(duration, t + 0.01);
			}
		}
		// DEBUG: std::ofstream("audio.raw", std::ios::binary).write(reinterpret_cast<char*>(&data[0]), data.size() * sizeof(float));
		// Filter the analyzer output data into QPainterPaths.
        m_allffts = analyzers[0].m_allffts;

        m_rate = analyzers[0].m_rate;
        m_fft_step = analyzers[0].processStep();
        // We need to create an image corresponding to the fft data
        unsigned fftx = analyzers[0].m_allffts.size();
        unsigned ffty = FFT_VIZ;
        FFTImage = QImage(fftx,ffty, QImage::Format_RGB32);
        FFTImage.fill(QColor(0,0,255));
        x=0;
        qDebug() <<"m_rate"<<m_rate<<"mpeg.audioqueue.rate"<<mpeg.audioQueue.getRate();
        qDebug() << "VizFFT info m_rate"<<m_rate<< "num_fft_samples" << analyzers[0].m_allffts.size() ;
        qDebug() <<"Calced fft samples:" << analyzers[0].m_allffts.size()*512 << "Duration = "<< duration << "Calcduration = " <<analyzers[0].m_allffts.size()*512.0/m_rate;

        //std::list<std::vector<double>>::iterator it_hf = analyzers[0].m_allffts_hf.begin();
        //for (std::list<std::vector<double>>::iterator itfft = analyzers[0]., itfftend = analyzers[0].m_allffts.end(); itfft != itfftend; itfft++){
        int first = 1;
        for (std::list<std::vector<double>>::iterator itfft = analyzers[0].m_allffts.begin(), itfftend = analyzers[0].m_allffts.end(); itfft != itfftend; itfft++){


            for (unsigned y = 0; y < ffty ; y++){
                double dbvalue = 96.0+itfft->at(y); // will peak at like 90
                if (first){
                    qDebug()<< "pitchviz create image y"<<y<<"dbvalue"<<dbvalue;

                }
                float scaledb = (std::min(1.0,std::max(0.0, (dbvalue-10.0)/(50.0))));// good for fft
                //int scaledb = 255.0*(std::min(1.0,std::max(0.0, (dbvalue -100.0)/60.0))); // for eac
                int r = 0;
                int g = 0;
                int b = 0;
                if (scaledb<0.25){
                    r=0;
                    g=0;
                    b=255*4*scaledb;
                }else if (scaledb < 0.5){
                    r=0;
                    g = 255*(4*(scaledb-0.25));
                    b = 255*(4*(0.5-scaledb));

                }else if (scaledb < 0.75){
                    r = (scaledb -0.5)*4*255;
                    g = (0.75-scaledb)*4*255;
                    b = 0;

                }else{
                    r = 255;
                    g = (scaledb-0.75)*4*255;
                    b = (scaledb-0.75)*4*255;

                }

                //double dbvalue_hf = 96.0+it_hf->at(y/4);
                //int scaledb_hf = 255.0*(std::min(1.0,std::max(0.0, (dbvalue_hf -20.0)/40.0)));

                FFTImage.setPixelColor(x,ffty-y-1,QColor(r,g,b));
                //FFTImage.setPixelColor(x,ffty-y-1,QColor(scaledb,scaledb,0));

            }
            //std::advance(it_hf, 1);
            x++;
            first = 0;
        }

		std::vector<Analyzer::Moments::const_iterator> mit(channels), mend(channels);
		for (unsigned ch = 0; ch < channels; ++ch) {
			Analyzer::Moments const& moments = analyzers[ch].getMoments();
			mit[ch] = moments.begin();
			mend[ch] = moments.end();
		}
		while (mit[0] != mend[0]) {
			for (unsigned ch = 0; ch < channels; ++mit[ch++]) {
				Moment::Tones const& tones = mit[ch]->m_tones;  // Take tones then move forward the iterator
				for (Moment::Tones::const_iterator it2 = tones.begin(), it2end = tones.end(); it2 != it2end; ++it2) {
					if (it2->prev) continue;  // The tone doesn't begin at this moment, skip
					// Copy the linked list into vector for easier access and calculate max level
					std::vector<Tone const*> tones;
					for (Tone const* n = &*it2; n; n = n->next) { tones.push_back(n); }
					if (tones.size() < 3) continue;  // Too short tone, ignored
					PitchPath path(ch);
					double score = 0.0;
					Analyzer::Moments::const_iterator momit = mit[ch];
					// Store path used for rendering
					for (unsigned i = 0; i < tones.size(); ++i, ++momit) {
						float t = momit->m_time;
						float n = scale.getNote(tones[i]->freq);
						float level = level2dB(tones[i]->level);

                        score += tones[i]->level;
						path.fragments.push_back(PitchFragment(t, n, level));
					}
					QMutexLocker locker(&mutex);
					if (score > 1.0) paths.push_back(path);
				}
			}
		}
		analyzingSuccess = true;

	} catch (std::exception& e) {
		std::cerr << std::string("Error loading audio: ") + e.what() + '\n' << std::flush;
	}
	{
		QMutexLocker locker(&mutex);
		moreAvailable = true;
		position = duration;
	}

	// Start the renderer loop
	if (analyzingSuccess) renderer();
}

void PitchVis::paint(int x1, int y1, int x2, int y2)
{
	QMutexLocker locker(&mutex);
	m_x1 = x1; m_y1 = y1;
	m_x2 = x2; m_y2 = y2;

	// Wake the thread
	restart = true;
	condition.wakeOne();
}

void PitchVis::renderer() {
	forever {
		int x1, x2, y1, y2;
		{
			QMutexLocker locker(&mutex);
			if (quit) return;
			x1 = m_x1, x2 = m_x2, y1 = m_y1, y2 = m_y2;
		}

		// Rendering
		// QImage allows drawing in non-main/non-GUI thread
		QImage image(x2-x1, y2-y1, QImage::Format_ARGB32_Premultiplied);
        //image = FFTImage.scaled(x2-x1,y2-y1,Qt::IgnoreAspectRatio,Qt::FastTransformation);


        NoteGraphWidget *widget = qobject_cast<NoteGraphWidget*>(parent());
		if (!widget) continue;

        double indextotime = m_rate / m_fft_step; //approxx 86 ffts per sec
        unsigned start = (widget->px2s(x1))*indextotime; //this transforms
        unsigned end = widget->px2s(x2)*indextotime;
        QImage drawslice = FFTImage.copy(start,0,end-start,FFT_VIZ);
        image =drawslice.scaled(x2-x1,y2-y1,Qt::IgnoreAspectRatio,Qt::FastTransformation);


        QSettings settings; // Default QSettings parameters given in main()
		bool aa = settings.value("anti-aliasing", true).toBool();
        //x12 and y12 are the sizes of the "canvas"
		{
			QPainter painter(&image);
			painter.setCompositionMode(QPainter::CompositionMode_Source);
			if (aa) painter.setRenderHint(QPainter::Antialiasing);
			// Fill the background, otherwise the image will have all kinds of carbage
            //painter.fillRect(image.rect(), QColor(0,0,0,0));
            QPen pen;
            pen.setWidth(1);
			pen.setCapStyle(Qt::RoundCap);
            pen.setColor(QColor(0, 0, 255, 128));
            painter.setPen(pen);
            // ok now to render BPM small bars and large bars
            double graphstarttime = widget->px2s(x1);
            double graphendtime = widget->px2s(x2);
            int quarter_beat_index = 0;
            for (double beattime = widget->getDelay(); beattime < graphendtime; beattime += (60.0/widget->getBPM())/4.0 ){
                if (beattime > graphstarttime){
                    int beatxpos = widget->s2px(beattime) -x1 ;
                    if (quarter_beat_index % 4 == 0) {
                        pen.setColor(QColor(255, 255, 255, 255));
                        painter.setPen(pen);

                    }else {
                        pen.setColor(QColor(0, 0, 100, 128));
                        painter.setPen(pen);

                    }

                    painter.drawLine(beatxpos,y1,beatxpos,y2);

                    //qDebug() << "Beat line: " << beatxpos << " time"<<beattime;
                }
                quarter_beat_index ++;
            }

            qDebug() << "graphstarttime: " << graphstarttime << " graphendtime"<<graphendtime << "m_rate"<<m_rate<<"x1"<<x1<<"x2"<<x2<<"y1"<<y1<<"y2"<<y2;


			PitchVis::Paths const& paths = getPaths();
			for (PitchVis::Paths::const_iterator it = paths.begin(), itend = paths.end(); it != itend; ++it) {
                break; //disable path rendering as its ugly
				PitchPath::Fragments const& fragments = it->fragments;
				int oldx, oldy;
				// Only render paths in view
				if (widget->s2px(fragments.back().time) < x1) continue;
				else if (widget->s2px(fragments.front().time) > x2) break;
				// Iterate through the path points
				for (PitchPath::Fragments::const_iterator it2 = fragments.begin(), it2end = fragments.end(); it2 != it2end; ++it2) {
					// TODO: Take y-size into account (change also the paint calls in NoteGraphWidget)
					int x = widget->s2px(it2->time) - x1;
					int y = widget->n2px(it2->note);
					if (m_visId == 0)
						pen.setColor(QColor(32 + 64 * it->channel, clamp<int>(127 + it2->level, 32, 255), 32, 128));
					else
						pen.setColor(QColor(clamp<int>(127 + it2->level, 32, 255), 32, 32 + 32 * it->channel, 100));
					painter.setPen(pen);
                    //qDebug() << "Drawing note  "<< it2->note << " str "<<qUtf8Printable(scale.getNoteStr(it2->note)) <<" freq " << scale.getNoteFreq(it2->note) <<" level "<< it2->level << " at time "<< it2->time;
					if (it2 != fragments.begin()) painter.drawLine(oldx, oldy, x, y);
					oldx = x; oldy = y;
				}
			}
		}

		// Send the image
		// This is actually delivered by the reciever's event loop thread, and not called directly from here
		emit renderedImage(image, QPoint(x1, y1), m_visId);

		mutex.lock();
		// If nothing to do, sleep here
		if (!restart) condition.wait(&mutex);
		restart = false;
		mutex.unlock();
	}
}

int PitchVis::guessNote(double begin, double end, int note) {
	const unsigned scoreSz = 48;
	double score[scoreSz] = {};
	if (note >= 0 || note < 48) score[note] = 10.0;  // Slightly prefer the current note
	// Score against paths
	for (PitchVis::Paths::const_iterator it = paths.begin(), itend = paths.end(); it != itend; ++it) {
		PitchPath::Fragments const& fragments = it->fragments;
		// Discard paths completely outside the window
		if (fragments.back().time < begin) continue;
		if (fragments.front().time > end) break;
		for (PitchPath::Fragments::const_iterator it2 = fragments.begin(), it2end = fragments.end(); it2 != it2end; ++it2) {
			// Discard path points outside the window
			if (it2->time < begin) continue;
			if (it2->time > end) break;
			unsigned n = round(it2->note);
			if (n < scoreSz) score[n] += 100 + it2->level;
		}
	}
	// Return the idx with best score
	return std::max_element(score + 1, score + scoreSz) - score;
}

