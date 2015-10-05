/***************************************************************************
 * Authors:     Vahid Abrishami (vabrishami@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <vector>
#include <sstream>
#include <fstream>
#include <time.h>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/video/video.hpp"
#include "opencv2/gpu/gpu.hpp"

#include <data/multidim_array.h>
#include <data/xmipp_image.h>
#include <data/normalize.h>
#include <data/xmipp_fftw.h>


using namespace std;
//using namespace cv;
#ifdef GPU
using namespace cv::gpu;
#endif
class ProgOpticalAligment: public XmippProgram
{

public:
    FileName fname, foname, gianRefFilename, darkRefFilename;
    MultidimArray<double> gainImage, darkImage;
    int winSize, gpuDevice, fstFrame, lstFrame;
    int psdPieceSize;
    bool doAverage, psd, saveCorrMovie;
    bool gainImageCorr, darkImageCorr;

    void defineParams()
    {
        addUsageLine ("Align moviews using optical flow");
        addParamsLine("     -i <inMoviewFnName>          : input movie File Name");
        addParamsLine("     -o <outAverageMoviewFnName>  : output aligned micrograhp File Name");
        addParamsLine("     [--nst <int=0>]     : first frame used in alignment (0 = first frame in the movie");
        addParamsLine("     [--ned <int=0>]     : last frame used in alignment (0 = last frame in the movie ");
        addParamsLine("     [--winSize <int=150>]     : window size for optical flow algorithm");
        addParamsLine("     [--simpleAverage]: if we want to just compute the simple average");
        addParamsLine("     [--psd]             : save raw PSD and corrected PSD");
        addParamsLine("     [--ssc]             : save corrected stack");
        addParamsLine("     [--gain <gainReference>]             : gain reference");
        addParamsLine("     [--dark <darkReference>]             : dark reference");
#ifdef GPU

        addParamsLine("     [--gpu <int=0>]         : GPU device to be used");
#endif

    }
    void readParams()
    {
        fname     = getParam("-i");
        foname    = getParam("-o");
        if ((gainImageCorr = checkParam("--gain")))
        {
            gianRefFilename = getParam("--gain");
        }
        if ((darkImageCorr = checkParam("--dark")))
        {
            darkRefFilename = getParam("--dark");
        }
        fstFrame  = getIntParam("--nst");
        lstFrame  = getIntParam("--ned");
        winSize   = getIntParam("--winSize");
        doAverage = checkParam("--simpleAverage");
        psd = checkParam("--psd");
        saveCorrMovie = checkParam("--ssc");

#ifdef GPU

        gpuDevice = getIntParam("--gpu");
#endif

    }
    void run()
    {
        main2();
    }

    // Save a matrix which is generated by OpenCV
    int saveMat( const string& filename, const cv::Mat& M)
    {
        if (M.empty())
        {
            return 0;
        }
        ofstream out(filename.c_str(), ios::out|ios::binary);
        if (!out)
            return 0;

        int cols = M.cols;
        int rows = M.rows;
        int chan = M.channels();
        int eSiz = (M.dataend-M.datastart)/(cols*rows*chan);

        // Write header
        out.write((char*)&cols,sizeof(cols));
        out.write((char*)&rows,sizeof(rows));
        out.write((char*)&chan,sizeof(chan));
        out.write((char*)&eSiz,sizeof(eSiz));

        // Write data.
        if (M.isContinuous())
        {
            out.write((char *)M.data,cols*rows*chan*eSiz);
        }
        else
        {
            return 0;
        }
        out.close();
        return 1;
    }

    // Load a matrix which is generated by saveMat
    int readMat( const string& filename, cv::Mat& M)
    {
        ifstream in(filename.c_str(), ios::in|ios::binary);
        if (!in)
        {
            //M = NULL_MATRIX;
            return 0;
        }
        int cols;
        int rows;
        int chan;
        int eSiz;

        // Read header
        in.read((char*)&cols,sizeof(cols));
        in.read((char*)&rows,sizeof(rows));
        in.read((char*)&chan,sizeof(chan));
        in.read((char*)&eSiz,sizeof(eSiz));

        // Determine type of the matrix
        int type = 0;
        switch (eSiz)
        {
        case sizeof(char):
                        type = CV_8UC(chan);
            break;
        case sizeof(float):
                        type = CV_32FC(chan);
            break;
        case sizeof(double):
                        type = CV_64FC(chan);
            break;
        }

        // Alocate Matrix.
        M = cv::Mat(rows,cols,type,cv::Scalar(1));

        // Read data.
        if (M.isContinuous())
{
            in.read((char *)M.data,cols*rows*chan*eSiz);
        }
        else
        {
            return 0;
        }
        in.close();
        return 1;
    }

    // Converts a XMIPP MultidimArray to OpenCV matrix
    void xmipp2Opencv(const MultidimArray<double> &xmippArray, cv::Mat &opencvMat)
    {
        int h = YSIZE(xmippArray);
        int w = XSIZE(xmippArray);
        opencvMat = cv::Mat::zeros(h, w,CV_32FC1);
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
        opencvMat.at<float>(i,j) = DIRECT_A2D_ELEM(xmippArray,i,j);
    }

    // Converts an OpenCV matrix to XMIPP MultidimArray
    void opencv2Xmipp(const cv::Mat &opencvMat, MultidimArray<double> &xmippArray)
    {
        int h = opencvMat.rows;
        int w = opencvMat.cols;
        xmippArray.initZeros(h, w);
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
        DIRECT_A2D_ELEM(xmippArray,i,j) = opencvMat.at<float>(i,j);
    }

    // Converts an OpenCV float matrix to an OpenCV Uint8 matrix
    void convert2Uint8(cv::Mat opencvDoubleMat, cv::Mat &opencvUintMat)
    {
        double min,max;
        cv::minMaxLoc(opencvDoubleMat, &min, &max);
        opencvDoubleMat.convertTo(opencvUintMat, CV_8U, 255.0/(max - min), -min * 255.0/(max - min));
    }

    // Computes the average of a number of frames in movies
    void computeAvg(const FileName &movieFile, int begin, int end, cv::Mat &avgimg)
    {
        ImageGeneric movieStack;
        MultidimArray<double> imgNormal;
        int N=end-begin+1;

        movieStack.readMapped(movieFile,begin);
        movieStack().getImage(imgNormal);
        if (darkImageCorr)
            imgNormal-=darkImage;
        if (gainImageCorr)
            imgNormal/=gainImage;
        xmipp2Opencv(imgNormal, avgimg);
        for (int i=begin;i<end;i++)
        {
            movieStack.readMapped(movieFile,i+1);
            movieStack().getImage(imgNormal);
            if (darkImageCorr)
                imgNormal-=darkImage;
            if (gainImageCorr)
                imgNormal/=gainImage;
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imgNormal)
            avgimg.at<float>(i,j)+=DIRECT_A2D_ELEM(imgNormal,i,j);
        }
        avgimg/=double(N);
        imgNormal.clear();
        movieStack.clear();
    }
    void std_dev2(const cv::Mat planes[], const cv::Mat &flowx, const cv::Mat &flowy, Matrix1D<double> &meanStdDev)
    {
        double sumX=0, sumY=0 ;
        double sqSumX=0, sqSumY=0;
        double valSubtract;
        int h=flowx.rows;
        int w=flowx.cols;
        int n=h*w;
        for(int i=0;i<h;i++)
            for(int j=0;j<w;j++)
            {
                valSubtract=planes[0].at<float>(i,j)-flowx.at<float>(i,j);
                sumX+=valSubtract;
                sqSumX+=valSubtract*valSubtract;
                valSubtract=planes[1].at<float>(i,j)-flowy.at<float>(i,j);
                sumY+=valSubtract;
                sqSumY+=valSubtract*valSubtract;
            }
        meanStdDev(0)=sumX/double(n);
        meanStdDev(1)=sqrt(sqSumX/double(n)-meanStdDev(0)*meanStdDev(0));
        meanStdDev(2)=sumY/double(n);
        meanStdDev(3)=sqrt(sqSumY/double(n)-meanStdDev(2)*meanStdDev(2));
    }

    int main2()
    {

        MultidimArray<double> preImg, avgCurr, avgStep, mappedImg;
        MultidimArray<double> outputMovie;
        Matrix1D<double> meanStdev;
        ImageGeneric movieStack, movieStackNormalize;
        Image<double> II;
        MetaData MD; // To save plot information
        FileName motionInfFile, correctedPSDFile, rawPSDFile;
        ArrayDim aDim;

        // For measuring times (both for whole process and for each level of the pyramid)
        clock_t tStart, tStart2;

#ifdef GPU
        // Matrix that we required in GPU part
        GpuMat d_flowx, d_flowy, d_dest;
        GpuMat d_avgcurr, d_preimg, d_mapx, d_mapy;
#endif

        // Matrix required by Opencv
        cv::Mat flowx, flowy, mapx, mapy, flow, dest;
        cv::Mat flowxPre, flowyPre, flowxInBet, flowyInBet;// Using for computing the plot information
        cv::Mat avgcurr, avgstep, preimg, preimg8, avgcurr8;
        cv::Mat planes[]={flowx, flowy};
        cv::Scalar meanx, meany;
        cv::Scalar stddevx, stddevy;

        int imagenum, cnt = 2, div = 0;
        int h, w, idx, levelNum, levelCounter = 1;

        motionInfFile=foname.replaceExtension("xmd");
        std::string extension=fname.getExtension();
        if (extension=="mrc")
            fname+=":mrcs";
        movieStack.read(fname,HEADER);
        movieStack.getDimensions(aDim);
        imagenum = aDim.ndim;
        h = aDim.ydim;
        w = aDim.xdim;
        if (darkImageCorr)
        {
            II.read(darkRefFilename);
            darkImage=II();
        }
        if (gainImageCorr)
        {
            II.read(gianRefFilename);
            gainImage=II();
        }
        meanStdev.initZeros(4);
        avgcurr=cv::Mat::zeros(h, w,CV_32FC1);
        flowxPre=cv::Mat::zeros(h, w,CV_32FC1);
        flowyPre=cv::Mat::zeros(h, w,CV_32FC1);
#ifdef GPU

        // Object for optical flow
        FarnebackOpticalFlow d_calc;
        setDevice(gpuDevice);

        // Initialize the parameters for optical flow structure
        d_calc.numLevels=6;
        d_calc.pyrScale=0.5;
        d_calc.fastPyramids=true;
        d_calc.winSize=winSize;
        d_calc.numIters=1;
        d_calc.polyN=5;
        d_calc.polySigma=1.1;
        d_calc.flags=0;
#endif
        // Initialize variables with zero
        // Initialize the stack for the output movie
        if (saveCorrMovie)
            outputMovie.initZeros(imagenum, 1, h, w);
        tStart2=clock();
        // Compute the average of the whole stack
        fstFrame++; // Just to adapt to Li algorithm
        lstFrame++; // Just to adapt to Li algorithm
        psdPieceSize = 400; // Currently we set it as a constant
        if (lstFrame>=imagenum || lstFrame==1)
            lstFrame=imagenum;
        imagenum=lstFrame-fstFrame+1;
        levelNum=sqrt(double(imagenum));
        computeAvg(fname, fstFrame, lstFrame, avgcurr);
        // if the user want to save the PSD
        if (psd)
        {
            opencv2Xmipp(avgcurr, avgCurr);
            II() = avgCurr;
            II.write(foname);
            if (doAverage)
                rawPSDFile=foname.removeLastExtension()+"_corrected";
            else
                rawPSDFile=foname.removeLastExtension()+"_raw";
            std::cerr<<"The file name is"<<rawPSDFile<<std::endl;
            String args=formatString("--micrograph %s --oroot %s --dont_estimate_ctf --pieceDim %d --overlap 0.7",
                                     foname.c_str(), rawPSDFile.c_str(), psdPieceSize);
            String cmd=(String)" xmipp_ctf_estimate_from_micrograph "+args;
            std::cerr<<"Computing the raw FFT"<<std::endl;
            if (system(cmd.c_str())==-1)
                REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
            if (doAverage)
                return 0;
            else
                foname.deleteFile();
        }
        avgCurr.clear();
        cout<<"Frames "<<fstFrame<<" to "<<lstFrame<<" under processing ..."<<std::endl;

        while (div!=1)
        {
            div = int(imagenum/cnt);
            // avgStep to hold the sum of aligned frames of each group at each step
            avgstep=cv::Mat::zeros(h, w,CV_32FC1);

            cout<<"Level "<<levelCounter<<"/"<<levelNum<<" of the pyramid is under processing"<<std::endl;
            // Compute time for each level
            tStart = clock();
            idx = 0;

            // Check if we are in the final step
            if (div==1)
                cnt = imagenum;

            for (int i=0;i<cnt;i++)
            {
                //Just compute the average in the last step
                if (div==1)
                {
                    movieStack.readMapped(fname,i+1);
                    movieStack().getImage(preImg);
                    if (darkImageCorr)
                        preImg-=darkImage;
                    if (gainImageCorr)
                        preImg/=gainImage;
                    xmipp2Opencv(preImg, preimg);
                }
                else
                {
                    if (i==cnt-1)
                        computeAvg(fname, i*div+fstFrame, lstFrame, preimg);
                    else
                        computeAvg(fname, i*div+fstFrame, (i+1)*div+fstFrame-1, preimg);
                }
                //xmipp2Opencv(preImg, preimg);
                // Note: we should use the OpenCV conversion to use it in optical flow
                convert2Uint8(avgcurr,avgcurr8);
                convert2Uint8(preimg,preimg8);
#ifdef GPU

                d_avgcurr.upload(avgcurr8);
                d_preimg.upload(preimg8);
                d_calc(d_avgcurr, d_preimg, d_flowx, d_flowy);
                d_flowx.download(flowx);
                d_flowy.download(flowy);
                d_avgcurr.release();
                d_preimg.release();
                d_flowx.release();
                d_flowy.release();
#else

                calcOpticalFlowFarneback(avgcurr8, preimg8, flow, 0.5, 6, winSize, 1, 5, 1.1, 0);
                split(flow, planes);
#endif
                // Save the flows if we are in the last step
                if (div==1)
                {
                    if (i > 0)
                    {
                        std_dev2(planes,flowxPre,flowyPre,meanStdev);
                        size_t id=MD.addObject();
                        MD.setValue(MDL_OPTICALFLOW_MEANX, double(meanStdev(0)), id);
                        MD.setValue(MDL_OPTICALFLOW_MEANY, double(meanStdev(2)), id);
                        MD.setValue(MDL_OPTICALFLOW_STDX, double(meanStdev(1)), id);
                        MD.setValue(MDL_OPTICALFLOW_STDY, double(meanStdev(3)), id);
                        MD.write(motionInfFile, MD_APPEND);
                    }
                    flowxPre = planes[0].clone();
                    flowyPre = planes[1].clone();
                }
                for( int row = 0; row < planes[0].rows; row++ )
                    for( int col = 0; col < planes[0].cols; col++ )
                    {
                        planes[0].at<float>(row,col) += (float)col;
                        planes[1].at<float>(row,col) += (float)row;
                    }
#ifdef GPU

                {
                    d_mapx.upload(mapx);
                    d_mapy.upload(mapy);
                    d_preimg.upload(preimg);
                    remap(d_preimg,d_dest,d_mapx,d_mapy,cv::INTER_CUBIC);
                    d_dest.download(dest);
                    d_dest.release();
                    d_preimg.release();
                    d_mapx.release();
                    d_mapy.release();
                }
#else
                cv::remap(preimg, dest, planes[0], planes[1], cv::INTER_CUBIC);
#endif

                if (div==1 && saveCorrMovie)
                    mappedImg.aliasImageInStack(outputMovie, i);
                avgstep+=dest;
            }
            avgcurr=avgstep/cnt;
            cout<<"Processing level "<<levelCounter<<"/"<<levelNum<<" has been finished"<<std::endl;
            printf("Processing time: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
            cnt=cnt*2;
            levelCounter++;
        }
        opencv2Xmipp(avgcurr, avgCurr);
        II() = avgCurr;
        II.write(foname);
        printf("Total Processing time: %.2fs\n", (double)(clock() - tStart2)/CLOCKS_PER_SEC);
        if (saveCorrMovie)
        {
            II()=outputMovie;
            II.write(foname.replaceExtension("mrcs"));
        }

        if (psd)
        {
            Image<double> psdCorr, psdRaw;
            MultidimArray<double> psdCorrArr, psdRawArr;
            correctedPSDFile=foname.removeLastExtension()+"_corrected";
            String args=formatString("--micrograph %s --oroot %s --dont_estimate_ctf --pieceDim %d --overlap 0.7",
                                     foname.c_str(), correctedPSDFile.c_str(), psdPieceSize);
            String cmd=(String)" xmipp_ctf_estimate_from_micrograph "+args;
            std::cerr<<"Computing the corrected FFT"<<std::endl;
            if (system(cmd.c_str())==-1)
                REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
            psdRaw.read(rawPSDFile+".psd");
            psdCorr.read(correctedPSDFile+".psd");
            psdCorrArr=psdCorr();
            psdRawArr=psdRaw();
            for (size_t i=0;i<psdPieceSize;i++)
                for (size_t j=0;j<size_t(psdPieceSize/2);j++)
                    DIRECT_A2D_ELEM(psdCorrArr,i,j)=DIRECT_A2D_ELEM(psdRawArr,i,j);
            psdCorr()=psdCorrArr;
            psdCorr.write(correctedPSDFile+".psd");
            FileName auxFile=rawPSDFile.addExtension("psd");
            auxFile.deleteFile();
        }

        // Release the memory
        avgstep.release();
        preimg.release();
        avgcurr8.release();
        preimg8.release();
        flow.release();
        planes[0].release();
        planes[1].release();
        flowxPre.release();
        flowyPre.release();
        movieStack.clear();
        preImg.clear();
        avgCurr.clear();
        II.clear();
        return 0;
    }
};

int main(int argc, char *argv[])
{
    ProgOpticalAligment prm;
    prm.read(argc,argv);
    return prm.tryRun();
}
