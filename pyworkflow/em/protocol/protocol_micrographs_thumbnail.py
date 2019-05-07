# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.protocol.constants import STEPS_PARALLEL
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.em.data import Image

from .protocol_micrographs import ProtPreprocessMicrographs


class ProtMicrographThumbnail(ProtPreprocessMicrographs):
    """
    Simple protocol to compute the micrograph thumbnails
    for visualization and analysis purposes.

    The current implementation of this protocol requires eman2
    to be installed to generate each thumbnail.
    """
    _label = 'micrograph thumbnail'

    def __init__(self, *args, **kwargs):
        ProtPreprocessMicrographs.__init__(self, *args, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      important=True,
                      label='Input micrographs',
                      help='Select a set of micrographs to generate the '
                           'thumbnails.')

        form.addParam('scaleFactor', params.IntParam, default=6,
                      label='Scale factor',
                      help='How many times (integer factor) do you want to '
                           'scale the micrograph to create the thumbnail. ')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        inputMics = self.inputMicrographs.get()
        deps = [self._insertFunctionStep("createThumbnailStep",
                                         mic.getFileName(), mic.getObjId())
                for mic in inputMics]

        self._insertFunctionStep("createOutputStep", prerequisites=deps)

    # --------------------------- STEPS functions ---------------------------------------------------
    def createThumbnailStep(self, micFn, micId):
        outputFn = self._getOutputMicThumbnail(micId)
        args = "%s %s " % (micFn, outputFn)
        args += "--fouriershrink %s --process normalize" % self.scaleFactor.get()

        import eman2
        pwutils.runJob(self._log, eman2.Plugin.getProgram('e2proc2d.py'), args,
                       env=eman2.Plugin.getEnviron())

    def createOutputStep(self):
        inputMics = self.inputMicrographs.get()
        outputMics = self._createSetOfMicrographs()
        outputMics.copyInfo(inputMics)
        for mic in inputMics:
            # Set the thumbnail image
            mic.thumbnail = Image(self._getOutputMicThumbnail(mic.getObjId()))
            outputMics.append(mic)

        self._defineOutputs(outputMicrographs=outputMics)
        self._defineTransformRelation(self.inputMicrographs, outputMics)

    # --------------------------- UTILS functions ----------------------------

    def _getOutputMicThumbnail(self, micId):
        return self._getExtraPath('mic%06d_thumbnail.png' % micId)