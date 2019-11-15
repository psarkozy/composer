#include "pitch.hh"

#include "libda/fft.hpp"
#include <cmath>
#include <numeric>

#include "fft_new.hpp"
#include <QDebug>

static const unsigned FFT_P = 14; //14 is best for me // FFT size setting, will use 2^FFT_P sample FFT

static const unsigned FFT_P_HF = 10;  // FFT si
static const std::size_t FFT_N_HF = 1 << FFT_P_HF;  // FFT size in samples


static const unsigned ZERO_PADDING_FACTOR = 1;
static const std::size_t FFT_N = 1 << FFT_P;  // FFT size in samples
static const std::size_t FFT_STEP = 512;  // Step size in samples, should be <= 0.25 * FFT_N. Low values cause high CPU usage.
static const double FFT_VIZ_MINFREQ = 67.0;
static const double FFT_VIZ_MAXFREQ = 989.0;

// Limit the range to avoid noise and useless computation
static const double FFT_MINFREQ = 45.0;
static const double FFT_MAXFREQ = 3000.0;

Tone::Tone(): freq(), level(), prev(), next() {
	for (std::size_t i = 0; i < MAXHARM; ++i) harmonics[i] = 0.0;
}

bool Tone::operator==(double f) const {
	return std::abs(freq / f - 1.0) < 0.06;  // Half semitone
}

Analyzer::Analyzer(double rate, std::string id):
  m_rate(rate),
  m_id(id),
  m_window(FFT_N),
  m_window_hf(FFT_N_HF),
  m_fftLastPhase(FFT_N / 2),
  m_oldfreq(0.0)
{
  	// Hamming window
	for (size_t i=0; i < FFT_N; i++) {
		m_window[i] = 0.53836 - 0.46164 * std::cos(2.0 * M_PI * i / (FFT_N - 1));
	}
    for (int i = 0; i < FFT_N; i++) { //Hann Window
        if (i < FFT_N/ZERO_PADDING_FACTOR)
        m_window[i] = 0.5 * (1 - cos(2*3.141592*i/(FFT_N/ZERO_PADDING_FACTOR)));
        else
        m_window[i] = 0.0;
    }
    for (int i = 0; i < FFT_N_HF; i++) { //Hann Window
        m_window_hf[i] = 0.5 * (1 - cos(2*3.141592*i/(FFT_N_HF)));
    }
}

unsigned Analyzer::processSize() const { return FFT_N; }
unsigned Analyzer::processStep() const { return FFT_STEP; }

void Analyzer::calcFFT(float* pcm) {

	m_fft = da::fft<FFT_P>(pcm, m_window);
    m_fft_hf = da::fft<FFT_P_HF>(pcm,m_window_hf);
    /*Complex cdata[FFT_N];

    for (int i = 0; i < m_window.size();i++){
        cdata[i].real(pcm[i]*m_window[i]);
        cdata[i].imag(0.0);
    }
    CArray fftdata(cdata, FFT_N); //sizeof complex

    fft2(fftdata);

    //Fourier fdata;
    m_fft.resize(m_window.size());
    for (int i = 0; i < m_window.size();i++){
        Complex compo = fftdata[i];
        m_fft[i].imag(fftdata[i].imag());
        m_fft[i].real(fftdata[i].real());
    }
    //m_fft = fftdata;
*/

}

namespace {
	bool sqrLT(float a, float b) { return a * a < b * b; }
	bool matchFreq(double f1, double f2) {
		return std::abs(f1 / f2 - 1.0) < 0.06;
	}
}

void Combo::combine(Peak const& p) {
	freq += p.level * p.freq;  // Multiplication for weighted average
	level += p.level;
}

bool Combo::match(double freqOther) const { return matchFreq(freq, freqOther); }

void Analyzer::calcTones() { // this is run once for each FFT_STEP
    // Precalculated constantsv

    const double freqPerBin = m_rate / (FFT_N/2) ; //BECAUSE Maxfreq is FFT_N/r
    //IMPORTANT: freqeuncies CHANGE with mono/stereo tracks!
    const double freqPerBin_hf = m_rate / (FFT_N_HF/2)  ; //BECAUSE Maxfreq is FFT_N/r
	const double phaseStep = 2.0 * M_PI * FFT_STEP / FFT_N;
	const double normCoeff = 1.0 / FFT_N;
	// Limit frequency range of processing
	const size_t kMin = std::max(size_t(3), size_t(FFT_MINFREQ / freqPerBin));
	const size_t kMax = std::min(FFT_N / 2, size_t(FFT_MAXFREQ / freqPerBin));
    m_peaks.resize(kMax);



    //qDebug() << "Analyzing";

    //create a 512 vector of floats for the fft data

	// Process FFT into peaks
    int highestlevel_bin = 0;
    float highestlevel = 0.0;

	for (size_t k = 1; k < kMax; ++k) {
		double level = normCoeff * std::abs(m_fft[k]);
        if (level > highestlevel){
            highestlevel = level;
            highestlevel_bin = k;
        }
		double phase = std::arg(m_fft[k]);
		// Use the reassignment method for calculating precise frequencies
		double delta = phase - m_fftLastPhase[k];
		m_fftLastPhase[k] = phase;
		delta -= k * phaseStep;  // Subtract the expected phase difference
		delta = remainder(delta, 2.0 * M_PI);  // Map the delta phase into +/- M_PI interval
		delta /= phaseStep;  // Calculate how much difference that makes during a step
		m_peaks[k].freqFFT = k * freqPerBin;  // Calculate the simple FFT frequency
		m_peaks[k].freq = (k + delta) * freqPerBin;  // Calculate the true frequency
		m_peaks[k].level = level;

	}
    if (highestlevel > 0.00001) {
        //qDebug() << "Highest level " << highestlevel << "frequency_bin"<<highestlevel_bin <<"freqency"<< highestlevel_bin *freqPerBin;
    }

    std::vector<double> visfftvec;
    visfftvec.resize(FFT_VIZ);
    std::vector<double> visfftvec_hf;
    visfftvec_hf.resize(FFT_VIZ);

    for (size_t n = 0; n < FFT_VIZ;n++){

        size_t bin_id = (int) floor((FFT_VIZ_MINFREQ/freqPerBin)*pow(2.0,n/(FFT_VIZ/4.0)));

        visfftvec[n] = level2dB(normCoeff * std::abs(m_fft[bin_id]));
        bin_id = (int) floor((FFT_VIZ_MAXFREQ/freqPerBin_hf)*pow(2.0,n/(FFT_VIZ/4.0)));
        visfftvec_hf[n] = level2dB(normCoeff * std::abs(m_fft_hf[bin_id]));

        //qDebug() << " N " <<  n << " bin_id: "<<bin_id << "fft val="<< visfftvec[n]<< "known frequency" << freqPerBin*bin_id;
    }
    m_allffts.push_back(visfftvec);
    m_allffts_hf.push_back(visfftvec_hf);

    // Filter peaks and combine adjacent peaks pointing at the same frequency into one
	typedef std::vector<Combo> Combos;
	Combos combos;
	for (size_t k = kMin; k < kMax; ++k) {
		Peak const& p = m_peaks[k];
		bool ok = p.level > 1e-3 && p.freq >= FFT_MINFREQ && p.freq <= FFT_MAXFREQ && std::abs(p.freqFFT - p.freq) < freqPerBin;
		if (!ok) continue;
		// Do we need to add a new Combo (rather than using the last one)?
		if (combos.empty() || !combos.back().match(p.freq)) combos.push_back(Combo());
		combos.back().combine(p);
	}
	// Convert sum frequencies into averages
	for (Combos::iterator it = combos.begin(), itend = combos.end(); it != itend; ++it) {
		it->freq /= it->level;
	}
	// Only keep a reasonable amount of strongest combos
	std::sort(combos.begin(), combos.end(), Combo::cmpByLevel);
	if (combos.size() > 30) combos.resize(30);
	// The order may not be strictly correct, fix it...
	std::sort(combos.begin(), combos.end(), Combo::cmpByFreq);
	// Try to combine combos into tones (collections of harmonics)
	Tones tones;
	for (Combos::const_iterator it = combos.begin(), itend = combos.end(); it != itend; ++it) {
		for (int div = 1; div <= 3; ++div) {  // Missing fundamental processing
			Tone tone;
			int plausibleHarmonics = 0;
			double basefreq = it->freq / div;
			if (basefreq < FFT_MINFREQ) break;  // Do not try any lower frequencies
			for (Combos::const_iterator harm = it; harm != itend; ++harm) {
				double ratio = harm->freq / basefreq;
				unsigned n = round(ratio);
				if (n > Tone::MAXHARM) break; // No more harmonics can be found
				if (std::abs(ratio - n) > 0.03) continue; // Frequency doesn't match
				if (n == 0) throw std::logic_error("combos not correctly sorted");
				if (n % div != 0) ++plausibleHarmonics;
				double l = harm->level;
				tone.harmonics[n - 1] += l;
				tone.level += l;
				tone.freq += l * harm->freq / n;  // The sum of all harmonics' fundies (weighted by l)
			}
			if (div > 1 && plausibleHarmonics < 3) continue;  // Not a proper missing fundamental
			tone.freq /= tone.level;  // Average instead of sum
			tones.push_back(tone);
		}
	}
	// Clean harmonics misdetected as fundamental
	tones.sort();
	for (Tones::iterator it = tones.begin(); it != tones.end(); ++it) {
		Tones::iterator it2 = it;
		++it2;
		while (it2 != tones.end()) {
			double ratio = it2->freq / it->freq;
			double diff = std::abs(ratio - round(ratio));
			bool erase = false;
			if (diff < 0.02 && it2->level < 2.0 * it->level) erase = true;  // Precisely harmonic and not much stronger than fundamental
			// Perform the action
			if (erase) it2 = tones.erase(it2); else ++it2;
		}
	}
	temporalMerge(tones);
}

void Analyzer::temporalMerge(Tones& tones) {
	if (!m_moments.empty()) {
		Tones& old = m_moments.back().m_tones;
		Tones::iterator it = tones.begin();
		// Iterate over old tones
		for (Tones::iterator oldit = old.begin(); oldit != old.end(); ++oldit) {
			// Try to find a matching new tone
			while (it != tones.end() && *it < *oldit) ++it;
			// If match found
			if (it != tones.end() && *it == *oldit) {
				// Link together the old and the new tones
				oldit->next = &*it;
				it->prev = &*oldit;
			}
		}
	}
	
    m_moments.push_back(Moment(m_moments.size() * processStep() / m_rate));
	m_moments.back().stealTones(tones);  // No pointers are invalidated
}

Moment::Moment(double t): m_time(t) {}

void Moment::stealTones(Tones& tones) {
	m_tones.swap(tones);
}

