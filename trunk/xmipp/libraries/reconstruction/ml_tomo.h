/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.uam.es (2004)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef _MLTOMO_H
#define _MLTOMO_H

#include <data/fftw.h>
#include <data/fft.h>
#include <data/args.h>
#include <data/funcs.h>
#include <data/selfile.h>
#include <data/docfile.h>
#include <data/image.h>
#include <data/geometry.h>
#include <data/filters.h>
#include <data/mask.h>
#include "ctf.h"
#include "sampling.h"
#include "symmetries.h"
#include "symmetrize.h"
#include <pthread.h>
#include <vector>

#define SIGNIFICANT_WEIGHT_LOW 1e-8
#define SMALLANGLE 2.75
#define MLTOMODATALINELENGTH 10
class Prog_ml_tomo_prm;

// Thread declaration
void * threadMLTomoExpectationSingleImage( void * data );

// This structure is needed to pass parameters to threadMLTomoExpectationSingleImage
typedef struct{
    int thread_id;
    int thread_num;
    Prog_ml_tomo_prm *prm;
    SelFile *SF;
    int *iter;
    double *wsum_sigma_noise;
    double *wsum_sigma_offset;
    double *sumfracweight;
    double *LL;
    std::vector<Matrix3D<double> > *wsumimgs;
    std::vector<Matrix3D<double> > *wsumweds;
    std::vector<VolumeXmippT<double> > *Iref;
    std::vector<Matrix1D<double> > *docfiledata;
    Matrix1D<double> *sumw;
} structThreadExpectationSingleImage ;

/**@defgroup ml_tomo ml_align2d (Maximum likelihood in 2D)
   @ingroup ReconsLibraryPrograms */
//@{
/** ml_tomo parameters. */
class Prog_ml_tomo_prm
{
public:
    /** Filenames reference selfile/image, fraction docfile & output rootname */
    FileName fn_sel, fn_ref, fn_root, fn_frac, fn_sym, fn_missing, fn_doc;
    /** Command line */
    std::string cline;
    /** Sigma value for expected pixel noise */
    double sigma_noise;
    /** sigma-value for origin offsets */
    double sigma_offset;
    /** Vector containing estimated fraction for each model */
    Matrix1D<double> alpha_k;
    /** Flag whether to fix estimates for model fractions */
    bool fix_fractions;
    /** Flag whether to fix estimate for sigma of origin offset */
    bool fix_sigma_offset;
    /** Flag whether to fix estimate for sigma of noise */
    bool fix_sigma_noise;
    /** Starting iteration */
    int istart;
    /** Number of iterations to be performed */
    int Niter;
    /** dimension of the images */
    int oridim, dim, dim3, hdim;
    double ddim3;
    /** Number of reference images */
    int nr_ref;
    /** Keep angles from docfile in generation of random subset averages */
    bool do_keep_angles;
    /** Total number of experimental images */
    int nr_exp_images;
    /** Sum of squared amplitudes of the references */
    std::vector<double> A2, corrA2;
    /** Verbose level:
        1: gives progress bar (=default)
        0: gives no output to screen at all */
    int verb;
    /** Stopping criterium */
    double eps;
    /** SelFile images (working, test and reference set) */
    SelFile SF, SFr;
    /** Vector for images to hold references (new & old) */
    std::vector < VolumeXmippT<double> > Iref, Iold, Iwed;
    /** Matrices for calculating PDF of (in-plane) translations */
    Matrix3D<double> P_phi, Mr2;
    /** Flag for generation of initial models from random subsets */
    bool do_generate_refs;

    /** Optimal refno and angno from previous iteration */
    std::vector<int> imgs_optrefno, imgs_optangno;
    std::vector<double> imgs_trymindiff, miss_nr_pixels;
    /** Number for missing data structure group */
    std::vector<int> imgs_missno;
    /** For initial guess of mindiff */
    double trymindiff_factor;
    /** Local search angular distance */
    double ang_search;
    /** Perturb angular sampling */
    bool do_perturb;
    /** Low-pass filter at FSC=0.5 resolution in each step */
    bool do_filter;

    /// Internally store all missing wedges or re-compute on the fly?

    /** Use missing wedges, cones or pyramids */
    bool do_missing, do_wedge, do_pyramid, do_cone;
    /** Do imputaion-like algorithm? */
    bool do_impute;
    /** Do maximum-likelihood algorithm? */
    bool do_ml;
    /** Threshold for no-imputation algorithm */
    double noimp_threshold;
    /** Number of missing data structures */
    int nr_miss;
    /** Maximum resolution (dig.freq.) */
    double maxres, scale_factor;
    Matrix3D<double> fourier_mask, fourier_imask, real_mask, real_omask;

    /* Adjust power spectra of all tilt series to average power spectrum */
    bool do_adjust_spectra;
    /** Arrays with power spectra for each tilt series */
    double *spectra_series;

    // Missing data information
#define MISSING_WEDGE_Y 0
#define MISSING_WEDGE_X 1
#define MISSING_PYRAMID 2
#define MISSING_CONE 3
    struct missing_info
    {
        int type;
        double thy0;
        double thyF;
        double thx0;
        double thxF;
    };
    typedef std::vector<missing_info> All_missing_info;
    All_missing_info all_missing_info;

    // Angular samnpling information
    struct angle_info
    {
        double rot;
        double tilt;
        double psi;
        int direction;
        Matrix2D<double> A;
    };
    typedef std::vector<angle_info> All_angle_info;
    /** Angular sampling  */
    double angular_sampling, psi_sampling;
    /** Vector with all angle combinations */
    All_angle_info all_angle_info;
    /** Number of angle combinations */
    int nr_ang;
    /** Pixel size */
    double pixel_size;

    /** Regularization parameters */
    double reg0, regF, reg_current;
    int reg_steps;

    /** Switch off SMALL_ANGLE addition (for phantoms) */
    bool no_SMALLANGLE;

    // sampling object
    XmippSampling mysampling;
    // For user-provided tilt range
    double tilt_range0, tilt_rangeF;
    // Symmetry setup
    int symmetry, sym_order;

    /** Threads */
    int threads;

    /** FFTW objects */
    XmippFftw transformer;

    /** debug flag */
    int debug;

public:
    /// Read arguments from command line
    void read(int argc, char **argv, bool ML3D = false);

    /// Show
    void show(bool ML3D = false);

    /// Usage for ML mode
    void usage();

    /// Extended Usage
    void extendedUsage(bool ML3D = false);

    /// Setup lots of stuff
    void produceSideInfo();

    /// Adjust power spectra
    void preparePowerSpectraAdjustment();

    /// Apply power spectra Adjustment to a single fourier transform
    void applyPowerSpectraAdjustment(Matrix3D<std::complex<double> > &M, int missno);

    /// Generate initial references from random subset averages
    void generateInitialReferences();

    /** Read reference images in memory & set offset vectors
        (This produce_side_info is Selfile-dependent!) */
    void produceSideInfo2(int nr_vols = 1);

    /// Calculate Angular sampling
    void calculateAngularSampling(int iter);

    /// Calculate probability density distribution for in-plane transformations
    void calculatePdfTranslations();

    /// Get binary missing wedge (or pyramid) 
    void getMissingRegion(Matrix3D<double> &Mmeasured,
                          Matrix2D<double> A,
                          const int missno);

    void maskSphericalAverageOutside(Matrix3D<double> &Min);

    // Resize a volume, based on the max_resol 
    // if down_scale=true: go from oridim to dim
    // if down_scale=false: go from dim to oridim
    void reScaleVolume(Matrix3D<double> &Min, bool down_scale=true);

    void postProcessVolume(VolumeXmipp &Vin, double &resolution);

    /// Fill vector of matrices with all rotations of reference
    void precalculateA2(std::vector< VolumeXmippT<double> > &Iref);

    /// ML-integration over all hidden parameters
    void expectationSingleImage(Matrix3D<double> &Mimg, int imgno, int missno,
                                std::vector<VolumeXmippT<double> > &Iref,
                                std::vector<Matrix3D<double> > &wsumimgs,
                                std::vector<Matrix3D<double> > &wsumweds,
                                double &wsum_sigma_noise, double &wsum_sigma_offset,
                                Matrix1D<double> &sumw, double &LL, double &dLL, 
                                double &fracweight, double &sumfracweight, double &trymindiff,
                                int &opt_refno, int &opt_angno, Matrix1D<double> &opt_offsets);

    /// Maximum constrained correlation search over all hidden parameters
    void maxConstrainedCorrSingleImage(Matrix3D<double> &Mimg, int imgno, int missno,
                                       std::vector<VolumeXmippT<double> > &Iref,
                                       std::vector<Matrix3D<double> > &wsumimgs,
                                       std::vector<Matrix3D<double> > &wsumweds,
                                       Matrix1D<double> &sumw, double &maxCC, double &sumCC,
                                       int &opt_refno, int &opt_angno, Matrix1D<double> &opt_offsets);

    /// Integrate over all experimental images
    void expectation(SelFile &SF, std::vector< VolumeXmippT<double> > &Iref, int iter,
                     double &LL, double &sumfracweight, DocFile &DFo,
                     std::vector<Matrix3D<double> > &wsumimgs,
                     std::vector<Matrix3D<double> > &wsumweds,
                     double &wsum_sigma_noise, double &wsum_sigma_offset,
                     Matrix1D<double> &sumw);

    /// Update all model parameters
    void maximization(std::vector<Matrix3D<double> > &wsumimgs,
                      std::vector<Matrix3D<double> > &wsumweds,
                      double &wsum_sigma_noise, double &wsum_sigma_offset, 
                      Matrix1D<double> &sumw, double &sumfracweight, 
                      double &sumw_allrefs, std::vector<Matrix1D<double> > &fsc, int iter);

    /// Calculate resolution by FSC
    void calculateFsc(Matrix3D<double> &M1, Matrix3D<double> &M2,
                      Matrix3D<double> &W1, Matrix3D<double> &W2,
                      Matrix1D< double >& freq, Matrix1D< double >& fsc, 
                      double &resolution);

    /// Apply regularization
    bool regularize(int iter);

    /// check convergence
    bool checkConvergence(std::vector<double> &conv);

    /// Write out reference images, selfile and logfile
    void writeOutputFiles(const int iter, DocFile &DFo,
                          std::vector<Matrix3D<double> > &wsumweds,
                          double &sumw_allrefs, double &LL, double &avefracweight,
                          std::vector<double> &conv, std::vector<Matrix1D<double> > &fsc);

};
//@}
#endif
