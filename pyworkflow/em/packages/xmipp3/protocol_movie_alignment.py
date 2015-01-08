# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Vahid Abrishami (vabrisahmi@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
In this module are protocol base classes related to EM Micrographs
"""

from os.path import join

from pyworkflow.object import String
from pyworkflow.protocol.params import IntParam, StringParam, BooleanParam, LEVEL_EXPERT, LEVEL_ADVANCED, EnumParam
from pyworkflow.utils.path import moveFile
import pyworkflow.em as em
from pyworkflow.em.protocol import ProtProcessMovies
from pyworkflow.gui.plotter import Plotter
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np

# Alignment methods enum
AL_OPTICAL = 0
AL_DOSEFGPU = 1
AL_DOSEFGPUOPTICAL = 2
AL_AVERAGE = 3

PLOT_CART = 0
PLOT_POLAR = 1
PLOT_POLARCART = 2
PLOT_CHOICES = ['cartesian', 'polar', 'both']

class ProtMovieAlignment(ProtProcessMovies):
    """ Aligns movies, from direct detectors cameras, into micrographs.
    """
    _label = 'movie alignment'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)

        form.addParam('alignMethod', EnumParam, choices=['optical flow', 'dosefgpu',
                                                         'dosefgpu + optical flow', 'average'],
                      label="Alignment method", default=AL_OPTICAL,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Method to use for alignment of the movies')
        group = form.addGroup('Common parameters')
        line = group.addLine('Used in alignment',
                            help='First and last frames used in alignment.\n'
                                  'The first frame in the stack is *0*.' )
        line.addParam('alignFrame0', IntParam, default=0, label='First')
        line.addParam('alignFrameN', IntParam, default=0, label='Last',
                      help='If *0*, use maximum value')
        group = form.addGroup('GPU parameters',condition="alignMethod==%d or alignMethod==%d "
                                                         " or alignMethod==%d "
                                                         % (AL_OPTICAL, AL_DOSEFGPUOPTICAL, AL_DOSEFGPU))
        group.addParam('doGPU', BooleanParam, default=False,
                      label="Use GPU (vs CPU)",
                      condition="alignMethod==%d or alignMethod==%d" % (AL_OPTICAL, AL_DOSEFGPUOPTICAL),
                      help="Set to true if you want the GPU implementation of Optical Flow")
        group.addParam('GPUCore', IntParam, default=1,
                      label="Choose GPU core",
                      condition="doGPU  or alignMethod==%d or alignMethod==%d  " % (AL_DOSEFGPU, AL_DOSEFGPUOPTICAL),
                      help="GPU may have several cores. Set it to one if you do not know what we are talking about")
        group = form.addGroup('Optical Flow parameters',condition="alignMethod==%d or alignMethod==%d " % (AL_OPTICAL, AL_DOSEFGPUOPTICAL))
        group.addParam('winSize', IntParam, default=150,
                      label="Window size", expertLevel=LEVEL_EXPERT,
                      help="Window size (shifts are assumed to be constant within this window).")
        #---------------------------------- DosefGPU Params--------------------------------
        group = form.addGroup('DosefGPU parameters',condition="alignMethod==%d or alignMethod==%d " % (AL_DOSEFGPU, AL_DOSEFGPUOPTICAL))
        line = group.addLine('Used in final sum',
                             help='First and last frames used in alignment.\n'
                                  'The first frame in the stack is *0*.' )
        line.addParam('sumFrame0', IntParam, default=0, label='First')
        line.addParam('sumFrameN', IntParam, default=0, label='Last',
                      help='If *0*, use maximum value')
        line = group.addLine('Crop offsets (px)')
        line.addParam('cropOffsetX', IntParam, default=0, label='X')
        line.addParam('cropOffsetY', IntParam, default=0, label='Y')
        line = group.addLine('Crop dimensions (px)',
                      help='How many pixels to crop from offset\n'
                           'If equal to 0, use maximum size.')
        line.addParam('cropDimX', IntParam, default=0, label='X')
        line.addParam('cropDimY', IntParam, default=0, label='Y')
        group.addParam('binFactor', IntParam, default=1,
                       label='Binning factor',
                       help='1x or 2x. Bin stack before processing.')
        group.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="""
-bft       150               BFactor in pix^2.
-pbx       96                Box dimension for searching CC peak.
-fod       2                 Number of frame offset for frame comparision.
-nps       0                 Radius of noise peak.
-sub       0                 1: Save as sub-area corrected sum. 0: Not.
-srs       0                 1: Save uncorrected sum. 0: Not.
-ssc       0                 1: Save aligned stack. 0: Not.
-scc       0                 1: Save CC Map. 0: Not.
-slg       1                 1: Save Log. 0: Not.
-atm       1                 1: Align to middle frame. 0: Not.
-dsp       1                 1: Save quick results. 0: Not.
-fsc       0                 1: Calculate and log FSC. 0: Not.
                      """)
        form.addParallelSection(threads=1, mpi=1)


    #--------------------------- STEPS functions ---------------------------------------------------
    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        micSet = self._createSetOfMicrographs()
        micSet.copyInfo(inputMovies)
        # Also create a Set of Movies with the alignment parameters
        movieSet = self._createSetOfMovies()
        movieSet.copyInfo(inputMovies)
        alMethod = self.alignMethod.get()
        for movie in self.inputMovies.get():
            micName = self._getNameExt(movie.getFileName(), 'mrc')
            metadataName = self._getNameExt(movie.getFileName(), 'xmd')
            plotName = self._getNameExt(movie.getFileName(), 'png')

            # Parse the alignment parameters and store the log files
            alignedMovie = movie.clone()
            alignedMovie.alignMetaData = String(self._getExtraPath(metadataName))
            print String(self._getExtraPath(plotName))
            alignedMovie.plotPolar = self._getExtraPath(plotName)
            movieCreatePlot(PLOT_POLAR, alignedMovie, True)
            movieSet.append(alignedMovie)
            print 'append has been done'


            mic = em.Micrograph()
            # All micrograph are copied to the 'extra' folder after each step
            mic.setFileName(self._getExtraPath(micName))
            mic.plotPolar = em.Image()
            mic.plotPolar.setFileName(self._getExtraPath(plotName))
            micSet.append(mic)


            # TODO: Methods for dosefgpu should be transferred to here
            """
            if alMethod == AL_DOSEFGPU:
                # Parse the alignment parameters and store the log files
                alignedMovie = movie.clone()
                logFile = self._getExtraPath(self._getLogFile(movie.getObjId()))
                import pyworkflow.em.packages.dosefgpu as dosefgpu
                alignment = dosefgpu.parseMovieAlignment(logFile)
                alignedMovie.setAlignment(alignment)
                movieSet.append(alignedMovie)
            """
        self._defineOutputs(outputMicrographs=micSet)
        self._defineOutputs(outputMovies=movieSet)
        self._defineSourceRelation(self.inputMovies.get(), micSet)
        """
        if alMethod == AL_DOSEFGPU:
            self._defineTransformRelation(inputMovies, micSet)
            self._defineOutputs(outputMovies=movieSet)
            self._defineTransformRelation(inputMovies, movieSet)
        """


    #--------------------------- UTILS functions ---------------------------------------------------

    def _processMovie(self, movieId, movieName, movieFolder):
        """ Process the movie actions, remember to:
        1) Generate all output files inside movieFolder (usually with cwd in runJob)
        2) Copy the important result files after processing (movieFolder will be deleted!!!)
        """

        program = self._getProgram()

        # Read the parameters
        #micName = self._getMicName(movieId)
        micName = self._getNameExt(movieName, 'mrc')
        print micName
        metadataName = self._getNameExt(movieName, 'xmd')
        firstFrame = self.alignFrame0.get()
        lastFrame = self.alignFrameN.get()
        gpuId = self.GPUCore.get() - 1
        alMethod = self.alignMethod.get()

        # For simple average execution
        if alMethod == AL_AVERAGE:
            command = '-i %(movieName)s -o %(micName)s' % locals()
            command += ' --nst %d --ned %d --simpleAverage' % (firstFrame, lastFrame)
            self.runJob(program, command, cwd=movieFolder)

        # For DosefGPU Execution (and combination with optical flow)
        if alMethod == AL_DOSEFGPU or alMethod == AL_DOSEFGPUOPTICAL:
            logFile = self._getLogFile(movieId)
            #gainFile = self.inputMovies.get().getGain()
            args = {'-crx': self.cropOffsetX.get(),
                    '-cry': self.cropOffsetY.get(),
                    '-cdx': self.cropDimX.get(),
                    '-cdy': self.cropDimY.get(),
                    '-bin': self.binFactor.get(),
                    '-nst': self.alignFrame0.get(),
                    '-ned': self.alignFrameN.get(),
                    '-nss': self.sumFrame0.get(),
                    '-nes': self.sumFrameN.get(),
                    '-gpu': gpuId,
                    '-flg': logFile,
                    }
            command = '%(movieName)s -fcs %(micName)s ' % locals()
            command += ' '.join(['%s %s' % (k, v) for k, v in args.iteritems()])
            if alMethod == AL_DOSEFGPUOPTICAL:
                program = 'dosefgpu_driftcorr'
                corrMovieName = self._getCorrMovieName(movieId)
                command += ' ' + '-fct %(corrMovieName)s -ssc 1 ' % locals()
            command += ' ' + self.extraParams.get()
            import pyworkflow.em.packages.dosefgpu as dosefgpu
            self.runJob(program, command, cwd=movieFolder,
                        env=dosefgpu.getEnviron())

        # For Optical Flow execution (and combination with DosefGPU)
        if alMethod == AL_OPTICAL or alMethod == AL_DOSEFGPUOPTICAL:
            winSize = self.winSize.get()
            if alMethod == AL_DOSEFGPUOPTICAL:
                program = 'xmipp_optical_alignment_gpu'
                corrMovieName = self._getCorrMovieName(movieId)
                command = '-i %(corrMovieName)s ' % locals()
                # Set to Zero for Optical Flow (output movie of dosefgpu)
                firstFrame = 0
                lastFrame = 0
            else:
                command = '-i %(movieName)s ' % locals()
            command += '-o %(micName)s --winSize %(winSize)d' % locals()
            command += ' --nst %d --ned %d' % (firstFrame, lastFrame)
            if self.doGPU:
                command += ' --gpu %d' % gpuId
            self.runJob(program, command, cwd=movieFolder)
            moveFile(join(movieFolder, metadataName), self._getExtraPath())

        # Move output micrograph and related information to 'extra' folder
        moveFile(join(movieFolder, micName), self._getExtraPath())
        if alMethod == AL_DOSEFGPU:
            # Copy the log file to have shifts information
            moveFile(join(movieFolder, logFile), self._getExtraPath())

    def _getProgram(self):
        alMethod = self.alignMethod.get()
        if alMethod == AL_AVERAGE:
            return 'xmipp_optical_alignment_cpu'
        if alMethod == AL_OPTICAL:
            if self.doGPU:
                return 'xmipp_optical_alignment_gpu'
            else:
                return 'xmipp_optical_alignment_cpu'
        elif alMethod == AL_DOSEFGPU:
            return 'dosefgpu_driftcorr'

  #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        numThreads = self.numberOfThreads;
        if numThreads>1:
            if self.doGPU:
                errors.append("GPU and Parallelization can not be used together")
        return errors

    def _citations(self):
        alMethod = self.alignMethod.get()
        if alMethod == AL_OPTICAL:
            return ['Abrishami2014a']
        if alMethod == AL_DOSEFGPU:
            return ['Li2013']
        if alMethod == AL_DOSEFGPUOPTICAL:
            return ['Abrishami2014a', 'Li2013']

    def _methods(self):
        methods = []
        alMethod = self.alignMethod.get()
        gpuId = self.GPUCore.get()
        if alMethod == AL_AVERAGE:
            methods.append('Aligning method: Simple average')
        if alMethod == AL_DOSEFGPU or alMethod == AL_DOSEFGPUOPTICAL:
            methods.append('Aligning method: DosefGPU')

        if alMethod == AL_OPTICAL or alMethod == AL_DOSEFGPUOPTICAL:
            methods.append('Aligning method: Optical Flow')
            methods.append('- Used a window size of: *%d*' % self.winSize.get())
            methods.append('- Used a pyramid size of: *6*')
            if self.doGPU:
                methods.append('- Used GPU *%d* for processing' % gpuId)

        return methods

    def _summary(self):
        firstFrame = self.alignFrame0.get()
        lastFrame = self.alignFrameN.get()
        summary = []
        summary.append('Number of input movies: *%d*' % self.inputMovies.get().getSize())
        if lastFrame == 0:
            summary.append('Frames used in alignment: *%d* to *%s*' % (firstFrame+1,'Last Frame'))
        else:
            summary.append('Frames used in alignment: *%d* to *%d*' % (firstFrame+1,lastFrame+1))

        return summary

def createPlots(plotType, protocol, movieId):

    movie = protocol.outputMovies[movieId]
    return movieCreatePlot(plotType, movie, False)

def movieCreatePlot(plotType, movie, saveFig):

    import xmipp
    meanX = []
    meanY = []
    stdX = []
    stdY = []
    colors = []
    gr = 255.0

    # if movie is None:
    #     return [self.errorMessage("Invalid movie id *%d*" % movieId,
    #                               title="Invalid input")]

    alignedMovie = movie.alignMetaData
    md = xmipp.MetaData(alignedMovie)
    colorDist = 255 / md.size()

    cartPosition = None
    polarPosition = None
    colorBarPosition = 122
    figureSize = (8, 6)

    if plotType == PLOT_CART:
        cartPosition = 121
    elif plotType == PLOT_POLAR:
        polarPosition = 121
    elif plotType == PLOT_POLARCART:
        cartPosition = 132
        polarPosition = 131
        colorBarPosition = 133

    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()

    # Plot the color bar
    ax = figure.add_subplot(colorBarPosition, aspect='equal', xlim=[0, 6])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colorBarX = np.array([1, 1])
    colorBarY = np.array([2, 4])

    # Read the shifts information from the Metadata
    for objId in md:
        meanX.append(md.getValue(xmipp.MDL_OPTICALFLOW_MEANX, objId))
        meanY.append(md.getValue(xmipp.MDL_OPTICALFLOW_MEANY, objId))
        stdX.append(md.getValue(xmipp.MDL_OPTICALFLOW_STDY, objId))
        stdY.append(md.getValue(xmipp.MDL_OPTICALFLOW_STDY, objId))
        colors.append((1, gr / 255.0, 0))
        ax.plot(colorBarX, colorBarY, c=(1, gr / 255.0, 0), linewidth=10)
        ax.text(2, np.mean(colorBarY), str(objId)+'-'+str(objId+1))
        colorBarY += 2
        gr -= colorDist
    area = (np.sqrt(np.power(np.asarray(stdX), 2) + np.power(np.asarray(stdY), 2)))*700

    # Plot in polar if needed
    if polarPosition:
        r = np.sqrt(np.power(np.asarray(meanX), 2) + np.power(np.asarray(meanY), 2))
        theta = np.arctan2(meanY, meanX) * 180 / np.pi
        ax = figure.add_subplot(polarPosition, projection='polar')
        ax.set_title('Polar representation')
        c = ax.scatter(theta, r, c=colors, s=area, cmap=plt.cm.hsv)
        c.set_alpha(0.75)
        ax.plot(theta, r, '-^')

    # Plot in cartesian if needed
    if cartPosition:
        ax = figure.add_subplot(cartPosition)
        ax.grid()
        ax.set_title('Cartesian representation')
        c = ax.scatter(np.asarray(meanX), np.asarray(meanY), c=colors, s=area, cmap=plt.cm.hsv)
        c.set_alpha(0.75)
        ax.plot(np.asarray(meanX), np.asarray(meanY), '-^')
    if saveFig:
        plotter.savefig(movie.plotPolar)
    return plotter