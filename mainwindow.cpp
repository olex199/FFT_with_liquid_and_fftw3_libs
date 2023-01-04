#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <iostream>

#include <complex.h>
#include <tgmath.h>

#include <liquid.h>
#include <fftw3.h>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    connect(ui->btnFFTExample, &QPushButton::clicked, this, &MainWindow::myFFt);
}

MainWindow::~MainWindow()
{
    delete ui;
}

#define NUM_POINTS 64

#define REAL 0
#define IMAG 1

void acquire_from_somewhere(fftw_complex* signal)
{
    /* Generate two sine waves of different frequencies and
     * amplitudes.
     */

    for (int i = 0; i < NUM_POINTS; ++i)
    {
        double theta = (double)i / (double)NUM_POINTS * M_PI;

        signal[i][REAL] = 1.0 * cos(10.0 * theta) +
                          0.5 * cos(25.0 * theta);

        signal[i][IMAG] = 1.0 * sin(10.0 * theta) +
                          0.5 * sin(25.0 * theta);
    }
}

void do_something_with(fftw_complex* result)
{
    for (int i = 0; i < NUM_POINTS; ++i)
    {
        double mag = sqrt(result[i][REAL] * result[i][REAL] +
                          result[i][IMAG] * result[i][IMAG]);

        //std::cout << "\n mag =  "<< mag;
        qDebug() << QString::number(i) << ": mag =  "<< mag;
    }
}

void MainWindow::myFFt()
{
    // liquid library **********************************************

    //unsigned int nfft = NUM_POINTS; // transform size
    int method = 0; // fft method (ignored)

    liquid_float_complex * xl = new liquid_float_complex[NUM_POINTS]{0};
    liquid_float_complex * yl = new liquid_float_complex[NUM_POINTS]{0};
    liquid_float_complex * zl = new liquid_float_complex[NUM_POINTS]{0};

    // initialize input
    for (uint i = 0; i < NUM_POINTS; i++)
    {
        double theta = (double)i / (double)NUM_POINTS * M_PI;

        xl[i].real(1.0 * cos(10.0 * theta) +
                  0.5 * cos(25.0 * theta));
        xl[i].imag(1.0 * sin(10.0 * theta) +
                  0.5 * sin(25.0 * theta));
    }

    // create fft plans
    fftplan pf = fft_create_plan(NUM_POINTS, xl, yl, LIQUID_FFT_FORWARD,  method);
    fftplan pr = fft_create_plan(NUM_POINTS, yl, zl, LIQUID_FFT_BACKWARD, method);

    // execute fft plans
    fft_execute(pf);
    fft_execute(pr);

    // destroy fft plans
    fft_destroy_plan(pf);
    fft_destroy_plan(pr);


//    qDebug() << "original signal, x[n]:\n";
//    for (uint i = 0; i < NUM_POINTS; i++)
//    {
//        //qDebug() << i << ". " << real(x[i]) << ", " << imag(x[i]);
//        qDebug() << i << ". " << x[i].real() << ", " << x[i].imag();
//    }

//    qDebug() << "\ny[n] = fft( x[n] ):\n";
//    for (uint i = 0; i < NUM_POINTS; i++)
//    {
//        //qDebug() << i << ". " << real(y[i]) << ", " << imag(y[i]);
//        qDebug() << i << ". " << y[i].real() << ", " << y[i].imag();
//    }

//    qDebug() << "\nz[n] = ifft( y[n] ):\n";
//    for (uint i = 0; i < NUM_POINTS; i++)
//    {
//        //qDebug() << i << ". " << real(z[i]) << ", " << imag(z[i]);
//        qDebug() << i << ". " << z[i].real()/(float)NUM_POINTS <<
//                    ", " << z[i].imag()/(float)NUM_POINTS;
//    }

    qDebug() << "done. __________________________________ \n";

    // fftw library **********************************************

    fftw_complex xf[NUM_POINTS];
    fftw_complex yf[NUM_POINTS];
    fftw_complex zf[NUM_POINTS];

    fftw_plan plan = fftw_plan_dft_1d(NUM_POINTS,
                                      xf,
                                      yf,
                                      FFTW_FORWARD,
                                      FFTW_ESTIMATE);

    fftw_plan back_plan = fftw_plan_dft_1d(NUM_POINTS,
                                      yf,
                                      zf,
                                      FFTW_BACKWARD,
                                      FFTW_ESTIMATE);

    for (int i = 0; i < NUM_POINTS; ++i)
    {
        double theta = (double)i / (double)NUM_POINTS * M_PI;

        xf[i][REAL] = 1.0 * cos(10.0 * theta) +
                      0.5 * cos(25.0 * theta);

        xf[i][IMAG] = 1.0 * sin(10.0 * theta) +
                      0.5 * sin(25.0 * theta);
    }
    //acquire_from_somewhere(xf);

    fftw_execute(plan);
    fftw_execute(back_plan);

    qDebug() << "\nliquid library: original signal, xl[n]:\t" <<
                "fftw library: original signal, xf[n]:\n";
    for (int i = 0; i < NUM_POINTS; ++i)
    {
        qDebug() << i << ". \t" <<
                    xl[i].real() << ", " << xl[i].imag() << "; \t" <<
                    (double)xf[i][REAL] << ", " << (double)xf[i][IMAG] <<
                    "; \t\t diff (must be 0): " <<
                    (double)(xl[i].real() - xf[i][REAL]) << ", " <<
                    (double)(xl[i].imag() - xf[i][IMAG]);
    }

    qDebug() << "\nliquid library: result signal, yl[n]:\t" <<
                "fftw library: result signal, yf[n]:\n";
    for (int i = 0; i < NUM_POINTS; ++i)
    {
        qDebug() << i << ". \t" <<
                    yl[i].real() << ", " << yl[i].imag() << "; \t" <<
                    (double)yf[i][REAL] << ", " << (double)yf[i][IMAG] <<
                    "; \t\t diff between liquid and fftw: " <<
                    yl[i].real() - (double)yf[i][REAL] << ", " <<
                    yl[i].imag() - (double)yf[i][IMAG];
    }

    qDebug() << "\nliquid library: back_result signal, zl[n]:\t" <<
                "fftw library: back_result signal, zf[n]:\n";
    for (int i = 0; i < NUM_POINTS; ++i)
    {
        qDebug() << i << ". \t" <<
                    zl[i].real()/(float)NUM_POINTS << ", " <<
                    zl[i].imag()/(float)NUM_POINTS << "; \t" <<
                    (double)zf[i][REAL]/(double)NUM_POINTS << ", " <<
                    (double)zf[i][IMAG]/(double)NUM_POINTS <<
                    "; \t\t diff between back and original: liquid:" <<
                    xl[i].real() - zl[i].real()/(float)NUM_POINTS << ", " <<
                    xl[i].imag() - zl[i].imag()/(float)NUM_POINTS << "; \tfftw:" <<
                    (double)(xf[i][REAL] - zf[i][REAL]/(double)NUM_POINTS) << ", " <<
                    (double)(xf[i][IMAG] - zf[i][IMAG]/(double)NUM_POINTS);
    }

//    do_something_with(signal);
//    do_something_with(result);

    delete [] xl;
    delete [] yl;
    delete [] zl;

    fftw_destroy_plan(plan);
    fftw_destroy_plan(back_plan);
}
