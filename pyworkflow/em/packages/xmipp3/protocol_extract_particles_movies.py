# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# *
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
This sub-package contains the XmippProtExtractParticles protocol
"""
 
import os

import pyworkflow.em.metadata as md
from pyworkflow.em.packages.xmipp3.convert import readSetOfParticles, xmippToLocation

from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.protocol import ProtExtractMovieParticles
# from pyworkflow.em.packages import xmipp3
# from pyworkflow.object import String, Float
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, IntParam, RelationParam,
                                        BooleanParam, Positive, FloatParam,
                                        EnumParam)
# from pyworkflow.em.protocol import  ProtExtractParticles
# from pyworkflow.em.data import SetOfParticles
from pyworkflow.em.constants import RELATION_CTF
from pyworkflow.utils.path import cleanPath, makePath, moveFile#removeBaseExt, replaceBaseExt
from pyworkflow.em.packages.xmipp3 import coordinateToRow, XmippMdRow
# 
# from convert import writeSetOfCoordinates, readSetOfParticles, micrographToCTFParam
# from xmipp3 import XmippProtocol
 
 
class XmippProtExtractMovieParticles(ProtExtractMovieParticles):
    """ Extract a set of Particles from each frame of a set of Movies.
    """
    _label = 'movie extract particles' 
     
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtExtractMovieParticles._defineParams(self, form)
        form.addParam('inputCoordinates', PointerParam, pointerClass='SetOfCoordinates',
                      important=True,
                      label='Input coordinates')
#         form.addParam('ctfRelations', RelationParam, allowsNull=True,
#                       relationName=RELATION_CTF, attributeName='getInputMovies',
#                       label='CTF estimation',
#                       help='Choose some CTF estimation related to input micrographs. \n'
#                            'CTF estimation is need if you want to do phase flipping or \n'
#                            'you want to associate CTF information to the particles.')
        form.addParam('boxSize', IntParam, default=0,
                      label='Particle box size', validators=[Positive],
                      help='In pixels. The box size is the size of the boxed particles, '
                      'actual particles may be smaller than this.')
        form.addParam('applyAlignment', BooleanParam, default=False,
                      label='Apply alignments to extract?')

        line = form.addLine('Frames range',
                      help='')
        line.addParam('firstFrame', IntParam, default=0,
                      label='First')
        line.addParam('lastFrame',  IntParam, default=0,
                      label='Last')
        
        form.addSection(label='Preprocess')
        form.addParam('doRemoveDust', BooleanParam, default=True, important=True,
                      label='Dust removal (Recommended)', 
                      help='Sets pixels with unusually large values to random values from a Gaussian '
                      'with zero-mean and unity-standard deviation.')
        form.addParam('thresholdDust', FloatParam, default=3.5, 
                      condition='doRemoveDust', expertLevel=LEVEL_ADVANCED,
                      label='Threshold for dust removal',
                      help='Pixels with a signal higher or lower than this value times the standard '
                      'deviation of the image will be affected. For cryo, 3.5 is a good value.'
                      'For high-contrast negative stain, the signal itself may be affected so '
                      'that a higher value may be preferable.')
        form.addParam('doInvert', BooleanParam, default=False,
                      label='Invert contrast', 
                      help='Invert the contrast if your particles are black over a white background.')
        form.addParam('doNormalize', BooleanParam, default=True,
                      label='Normalize (Recommended)', 
                      help='It subtract a ramp in the gray values and normalizes so that in the '
                      'background there is 0 mean and standard deviation 1.')
        form.addParam('normType', EnumParam, choices=['OldXmipp','NewXmipp','Ramp'], default=2, 
                      condition='doNormalize', expertLevel=LEVEL_ADVANCED,
                      display=EnumParam.DISPLAY_COMBO,
                      label='Normalization type', 
                      help='OldXmipp (mean(Image)=0, stddev(Image)=1).  \n  '
                           'NewXmipp (mean(background)=0, stddev(background)=1)  \n  '
                           'Ramp (subtract background+NewXmipp).  \n  ')
        form.addParam('backRadius', IntParam, default=-1, condition='doNormalize',
                      label='Background radius',
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                      'is set to 1. Radius for background circle definition (in pix.). '
                      'If this value is 0, then half the box size is used.', 
                      expertLevel=LEVEL_ADVANCED)

        form.addParallelSection(threads=0, mpi=0)
    
    #--------------------------- STEPS functions --------------------------------------------------
    def _processMovie(self, movieId, movieName, movieFolder):
        
        movieName = os.path.join(movieFolder, movieName)
        if self.doRemoveDust:
            self._runNoDust(movieName)
        
        boxSize = self.boxSize.get()
        # Read movie dimensions to iterate through each frame
        imgh = ImageHandler()
        x, y, z, n = imgh.getDimensions(movieName)
        
        first = self.firstFrame.get()
        last = self.lastFrame.get()
        
        if last == 0:
            last = n-1
        
        movie = self.inputMovies.get()[movieId]
        if self.applyAlignment and movie.hasAlignment():
            shifts = movie.getAlignment().getShifts()
        else:
            shifts = [0] * (2*n)
        
        stkIndex = 0
        movieStk = self._getMovieName(movieId, '.stk')
        movieMdFile = self._getMovieName(movieId, '.xmd')
        movieMd = md.MetaData()
        frameMd = md.MetaData()
        frameMdImages = md.MetaData()
        frameRow = md.Row()
         
        for frame in range(first, last+1):
            i = frame + 1
            # Get the frame shifts
            index = frame - first
            shiftX = shifts[2*index]
            shiftY = shifts[2*index+1]
            
            frameRoot = os.path.join(movieFolder, 'frame_%02d' % frame)
            frameName = frameRoot + '.mrc'
            frameMdFile = frameRoot + '.xmd'
            framePosFile = frameRoot + '_coordinates.xmd'
            coordinatesName = frameRoot + '_coordinates.xmd'
            
            hasCoordinates = self._writeXmippPosFile(movieId, movieName, coordinatesName, 
                                                     shiftX, shiftY)
            
            if hasCoordinates:
                self.info("Writing frame: %s" % frameName)
                imgh.convert(tuple([i, movieName]), frameName)
                 
                self.info("Extracting particles")
                frameImages = frameRoot + '_images'
                args = '-i %(frameName)s --pos %(coordinatesName)s -o %(frameRoot)s --Xdim %(boxSize)d' % locals()
                
                if self.doInvert:
                    args += " --invert"
                
                self.runJob('xmipp_micrograph_scissor', args)
                cleanPath(frameName)
                
                frameStk = frameRoot + '.stk'
                 
                self.info("Combining particles into one stack.")
                 
                frameMdImages.read(frameMdFile)
                frameMd.read('particles@%s' % framePosFile)
                frameMd.merge(frameMdImages)
                 
                for objId in frameMd:
                    stkIndex += 1
                    frameRow.readFromMd(frameMd, objId)
                    location = xmippToLocation(frameRow.getValue(md.MDL_IMAGE))
                    newLocation = (stkIndex, movieStk)
                    imgh.convert(location, newLocation)
                    
                    # Fix the name to be accesible from the Project directory
                    # so we know that the movie stack file will be moved
                    # to final particles folder
                    newImageName = '%d@%s' % newLocation
                    frameRow.setValue(md.MDL_IMAGE, newImageName)
                    frameRow.setValue(md.MDL_MICROGRAPH_ID, long(movieId))
                    frameRow.setValue(md.MDL_MICROGRAPH, str(movieId))
                    frameRow.setValue(md.MDL_FRAME_ID, long(frame))
                    frameRow.setValue(md.MDL_PARTICLE_ID, frameRow.getValue(md.MDL_ITEM_ID))
                    frameRow.writeToMd(movieMd, movieMd.addObject())
                movieMd.addItemId()
                movieMd.write(movieMdFile)
                cleanPath(frameStk)
        
        if self.doNormalize:
            self._runNormalize(movieStk)

    
    def _runNoDust(self, movieName):
        fnNoDust = movieName + "_noDust.stk"
        thresholdDust = self.thresholdDust.get() #TODO: remove this extra variable
        args=" -i %(movieName)s -o %(fnNoDust)s --bad_pixels outliers %(thresholdDust)f" % locals()
        self.runJob("xmipp_transform_filter", args)
        ImageHandler().convertStack(fnNoDust, movieName)
    
    def _runNormalize(self, stack):
        program = "xmipp_transform_normalize"
        args = "-i %(stack)s "
        normType = self.normType.get()
        bgRadius = self.backRadius.get()
        
        size, _, _ , _ = ImageHandler().getDimensions(stack)
        if bgRadius <= 0:
            bgRadius = int(size/2)
        
        if normType=="OldXmipp":
            args += "--method OldXmipp"
        elif normType=="NewXmipp":
            args += "--method NewXmipp --background circle %(bgRadius)d"
        else:
            args += "--method Ramp --background circle %(bgRadius)d"
        self.runJob(program, args % locals())
    
    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        particleSet = self._createSetOfMovieParticles()        
        particleSet.copyInfo(inputMovies)
          
        # Create a folder to store the resulting micrographs
         
#         particleFolder = self._getPath('particles')
#         makePath(particleFolder)
        mData = md.MetaData()
        mdAll = md.MetaData()
          
        for movie in inputMovies:
            movieId = movie.getObjId()
#             print "modieId: ", movieId
            movieName = self._getMovieName(movieId)
            movieStk = movieName.replace('.mrc', '.stk')
            movieMdFile = movieName.replace('.mrc', '.xmd')
             
            # Store particle stack and metadata files in final particles folder
            if os.path.exists(movieStk):
                mData.read(movieMdFile)
                mdAll.unionAll(mData)
              
        particleMd = self._getPath('movie_particles.xmd')
        mdAll.addItemId()
#         print "mdAll: ", mdAll
        mdAll.write(particleMd)
        readSetOfParticles(particleMd, particleSet, removeDisabled=False)
          
        self._defineOutputs(outputParticles=particleSet)
        self._defineSourceRelation(inputMovies, particleSet)
    
    #--------------------------- UTILS functions --------------------------------------------------
    def _writeXmippPosFile(self, movieId, movieName, coordinatesName, shiftX, shiftY):
        """ Create an Xmipp coordinates files to be extracted
        from the frames of the movie.
        """
        coordinates = self.inputCoordinates.get()
        micrograph = coordinates.getMicrographs()[movieId]
         
        if micrograph is None:
            raise Exception("Micrograph with id %d was not found in SetOfCoordinates!!!" % movieId)
         
        mData = md.MetaData()
        coordRow = XmippMdRow()
         
        for coord in coordinates.iterCoordinates(micrograph):
            coord.shiftX(1 * int(round(float(shiftX))))
            coord.shiftY(1 * int(round(float(shiftY))))
            coordinateToRow(coord, coordRow)
            coordRow.writeToMd(mData, mData.addObject())
             
        if mData.isEmpty():
            return False
        else:
            self.info("Writing coordinates metadata: %s, with shifts: %s %s" % (coordinatesName, shiftX, shiftY))
            mData.write('particles@' + coordinatesName)
            return True
    
    def _filterMovie(self, movieId, movieFn):
        """ Check if process or not this movie.
        In this case if there are not coordinates, do not process it.
        """
        coordinates = self.inputCoordinates.get()
        micrograph = coordinates.getMicrographs()[movieId]
         
        if micrograph is None:
            self.warning("Micrograph with id %d was not found in SetOfCoordinates!!!" % movieId)
            return False
         
        n = len([coord for coord in coordinates.iterCoordinates(micrograph)])
         
        if n == 0:
            self.warning("Skipping micrograph with id %d, not coordinates found!!!" % movieId)
            return False
         
        return True