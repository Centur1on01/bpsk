//#define Complex complex<double> // так не надо
#define _USE_MATH_DEFINES

#include <iostream>
#include <complex>
#include <string>
#include <fstream>
#include <cmath>
#include <deque>
#include <queue>
#include <random>
#include <vector>

using namespace std;

typedef complex<double> Complex;    // лучше так

// источник сигнала ФМ-2 (0 дБ)
class gen {
    complex<double> result;
    complex<double> temp;
    vector<complex<double>> a;
    ofstream text;
    ifstream textread;
    string pls = "";

public:
    gen()
    {
    }

    gen(string result) {
        a.resize(2);
        a[0].real(1);
        a[0].imag(0);
        a[1].real(-1);
        a[1].imag(0);
        text.open(result);
        textread.open(result);
        srand(time(0));
    }

    ~gen() {
        text.close();
        textread.close();
    }

    complex<double> next() {
        rand() % 2 == 0 ? temp = a[0] : temp = a[1];
        result = temp;
        temp.imag() < 0 ? pls = "" : pls = "+";
        text << temp.real() << pls << temp.imag() << "i" << endl;
        return result;
    }

    complex<double> getFromFile() {
        string line;
        if (getline(textread, line)) {
            double real, imag;
            char plus;
            char i;
            stringstream ss(line);
            ss >> real >> plus >> imag >> i;
            imag *= (plus == '+' ? 1 : -1);
            result.real(real);
            result.imag(imag);
        }
        return result;
    }

    void start() {
        //fseek(textread, 0L, SEEK_SET);
        textread.clear();
        textread.seekg(0, ios::beg);
    }

    bool eof() {
        return textread.eof();
    }
};
 



// !инициализацию необходимо выполнить, задавая требуемую мощность шума (в Вт или Дб)
class awgn {
    default_random_engine generator;
    normal_distribution<double> distr;

public:

    awgn(const double& m_db, const double& s_db)
    {
        //double m = pow(10, m_db / 20);
        //double s = pow(10, s_db / 20);
        distr = normal_distribution<double>(pow(10, m_db / 20), pow(10, s_db / 20));
        //distr = normal_distribution<double>(pow(10, m_db / 20), pow(10, s_db / 20));
    }

    double next() {
        return distr(generator);
    }

};

// Фильтр нижних частот
class LPF {
    int filterSize;
    double filterSum;
    double filteredSignal = 0;
    queue<double, deque<double>> q;
    deque<Complex> buf;

public:

    LPF(const int& filterSize = 8) : filterSize(filterSize), filterSum(0) {
        buf.resize(filterSize, {});
        for (int i = 0; i < filterSize; i++) {
            q.push(0);
        }
    }

    double filter(const double& signal) {
        filterSum -= q.front();
        q.pop();
        q.push(signal);
        filterSum += signal;
        filteredSignal = filterSum / filterSize;
        return filteredSignal;
    }

    double getFilteredSignal() {
        return filteredSignal;
    }

};

complex<double> jj(0.0, 1.0);

// Блок преобразования Фурье
class DFT {
    int N;
    Complex fourier;
    Complex signalTransformed = 0;

public:

    DFT(const int& N = 128) : N(N) {}

    // Возвращает сигнал после преобразования Фурье
    Complex transform(const vector<Complex>& signal, const int& i) {
        signalTransformed = 0;
        for (int j = 0; j < N; j++) {
            fourier = exp(-jj * 2.0 * M_PI * (double)i * (double)j / (double)N);
            signalTransformed += fourier * signal[j];
        }
        return signalTransformed;
    }

};

class decision {
    complex<double> result;

public:

    complex<double> makeDecision(complex<double> signal) {
        // y = -x ; imag = -real
        if (signal.imag() >= -signal.real()) {
            result.real(1);
            result.imag(0);
            return result;
        }
        else {
            result.real(-1);
            result.imag(0);
            return result;
        }
    }

    complex<double> makeReal(complex<double> signal) {
        if (signal.real() >= 0) {
            result.real(1);
            result.imag(0);
            return result;
        }
        else {
            result.real(-1);
            result.imag(0);
            return result;
        }
    }

};

int main()
{
    
    double snr_t = 0;
    ofstream output("theoretical_output.txt");

    // теоретичская кривая помехоустойчивости
    while (snr_t <= 15) { // цикл по значениям snr от 0 до 15 дБ
        double ber = 0.5 * erfc(sqrt(pow(10, snr_t / 10))); // вычисление вероятности ошибки бита
        output << snr_t << "\t" << ber << endl; // запись значения snr и ber в файл
        ++snr_t; // увеличение значения snr на 1 дБ
    }
    output.close(); 

    // !шум для заданного ОСШ
    double snr = 0;
    double signallevel = 0;
    double noiselevel = signallevel - snr;
    gen bitGenerator("result.txt");

    awgn noise(noiselevel, 0.0);
    decision decoder;
    ofstream outputp("outputp.txt"); 

    complex<double> jjj = { 1.0, 1.0 };
    // кривая помехоустойчивости для ФМ-2 с ОСШ от 0 до 15 дБ
    while (snr <= 15) { // цикл по значениям snr от 0 до 15 дБ
        int N = 100000; // количество передаваемых символов
        int errors = 0; // количество ошибок декодирования
        for (int i = 0; i < N; i++) {
            complex<double> signal = bitGenerator.next(); // генерация очередного символа
            complex<double> noisedSignal = signal + noise.next() * jjj; // добавление шума к сигналу
            complex<double> decodedSignal = decoder.makeDecision(noisedSignal); // декодирование сигнала
            if (decodedSignal != signal) ++errors; // подсчет ошибок декодирования
        }
        double ber = (double)errors / (double)N; // вычисление вероятности ошибки бита
        outputp << snr << "\t" << ber << endl; // запись значения snr и ber в файл
        cout << snr << '\t' << ber << endl;
        ++snr; // увеличение значения snr на 1 дБ
        //signallevel = signallevel + 1; // увеличение уровня сигнала на 1 дБ
        noiselevel = signallevel - snr; // вычисление уровня шума для нового значения snr
        awgn newNoise(noiselevel, 0.0); // обновление параметров объекта класса awgn
        noise = newNoise;
    }
    outputp.close();
    cout << "123" << endl;





    // старая реализация вычисления кривой помехоустойчивости и теоретических значений
    double Nbit = 1000;
    complex<double> tempbit;
    for (int i = 0; i < Nbit; i++) {
        tempbit = bitGenerator.next();
    }

    awgn noiseGenerator(noiselevel, 1.0);

    // кривая помехоустойчивости в matlab (0-15 дБ) - наложить теоретическую кривую для ФМ-2

    complex<double> sumResult;
    ofstream sumtext("sumresult.txt");
    string pls;

    decision decisionMaker;
    complex<double> signal;

    double EbN0_dB;
    vector<double> EbN0;
    vector<double> BER;
    double BERratio = 0;
    sumResult.imag(0);
    double BER_t;
    ofstream result_real("real.txt");
    double incorrect = 0, correct = 0;
    for (EbN0_dB = 0; EbN0_dB <= 15; EbN0_dB++) {
        incorrect = 0;
        correct = 0;
        EbN0.push_back(pow(10., EbN0_dB / 10.)); // pow(10.,EbN0_dB / 10)
        awgn noise(EbN0_dB, 0.0);
        //awgn noise(0, 1 / (sqrt(2) * EbN0[EbN0_dB]));
        bitGenerator.start();
        int it = 0;
        while (1) {
            if (it > Nbit - 1)
                break;
            signal = bitGenerator.getFromFile();
            sumResult.real(signal.real() + noise.next());
            if (signal == decisionMaker.makeReal(sumResult))
                ++correct;
            else ++incorrect;
            ++it;
            sumResult.imag() < 0 ? pls = "" : pls = "+";
            sumtext << sumResult.real() << pls << sumResult.imag() << "i" << endl;
        }
        BER_t = 0.5 * erfc(sqrt(EbN0[EbN0_dB]));
        BERratio = incorrect / (Nbit);
        cout << BERratio << "\t\t\t" << BER_t << endl;
        result_real << EbN0_dB << "\t" << BERratio << "\t\t\t" << BER_t << endl;
    }
    result_real.close();












    /*awgn signal(0.0, 1.0);
    LPF lpfilter(16);
    double f = 1000;
    double fs = 8000;
    double N = 1024;
    vector<complex<double>> signalvector;
    signalvector.resize(N);
    ofstream signalb4("signal.txt");
    for (int i = 0; i < N; ++i) {
        signalvector[i] = signal.next();
        signalb4 << i << "\t" << signalvector[i].real() << (signalvector[i].imag() >= 0 ? "+" : "") << signalvector[i].imag() << "i" << endl;
    }
    signalb4.close();
    DFT dft(N);

    vector<complex<double>> filteredsignalvector;
    filteredsignalvector.resize(N);
    ofstream signalf("signalf.txt");
    for (int i = 0; i < N; ++i) {
        filteredsignalvector[i] = lpfilter.filter(signalvector[i].real());
        signalf << i << "\t" << filteredsignalvector[i].real() << (filteredsignalvector[i].imag() >= 0 ? "+" : "") << filteredsignalvector[i].imag() << "i" << endl;
    }
    signalf.close();

    double nf = 0;
    double temp = 0;
    ofstream spectreOut("spectre.txt");
    ofstream spectreLPFOut("spectrelpf.txt");
    complex<double> spectre = 0;
    complex<double> spectreLPF = 0;
    for (int i = 0; i < N; i++) {
        if (i < N / 2) {
            temp = i + N / 2;
            spectre = dft.transform(signalvector, temp);
            spectreLPF = dft.transform(filteredsignalvector, temp);
        }
        else {
            temp = i - N / 2;
            spectre = dft.transform(signalvector, temp);
            spectreLPF = dft.transform(filteredsignalvector, temp);
        }

        nf = fs / N * (i - N / 2);
        spectreOut << nf << "\t" << abs(spectre) << endl;
        spectreLPFOut << nf << "\t" << abs(spectreLPF) << endl;

    }
    spectreOut.close();
    spectreLPFOut.close();*/

    std::system("pause");
}

