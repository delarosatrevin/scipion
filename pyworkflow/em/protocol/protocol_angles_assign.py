# **************************************************************************
# *
# * Authors:        Olivia Pfeil-Gardiner (zolivia@zedat.fu-berlin.de) [1]
#                   J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# *  [1] OPIC, Division of Structural Biology, University of Oxford
# *  [2] SciLifeLab, Stockholm University
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
from pyworkflow.em.protocol.protocol_particles import ProtProcessParticles
from pyworkflow.protocol.params import BooleanParam
from pyworkflow.em import metadata as md
import pyworkflow.em.packages.relion.convert
from random import randint

class ProtAnglesAssign(ProtProcessParticles):
    """Allows the user to assign a tilt angle of 90 deg and randomize the rotational angle of particles"""
    _label = 'assign angles'

    def __init__(self, **kwargs):
        ProtProcessParticles.__init__(self, **kwargs)

    #--------------- DEFINE param functions ---------------

    def _defineParams(self, form):

        ProtProcessParticles._defineParams(self, form)

    def _defineProcessParams(self,form):
        form.addParam('assignTilt90', BooleanParam, label='Assign out of plane angle (Tilt) to 90 deg?',
                      help='Assigns a value of 90 degrees to the tilt angle of each particle.')

        form.addParam('randomizeRot', BooleanParam,
                      label='Randomize the rotational angle around the particle axis (Rot)?',
                      help='Assigns random values to the rotational angles around the particle axis.')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions ---------------------------------------------------

    def createOutputStep(self):
        pass
#        particles = self.inputParticles.get()

      #  if self.assignTilt90:
     #       alignRow = md.Row()
    #        particles.setAlignment2D()
   #         for part in particles.iterItems():
  #              alignRow.setValue(md.RLN_ORIENT_TILT, 90)
 #               transform = rowToAlignment(alignRow, em.ALIGN_2D)
#                part.setTransform(transform)

        #if self.randomizeRot:
       #     alignRow = md.Row()
      #      particles.setAlignment2D()
     #       for part in particles.iterItems():
    #            alignRow.setValue(md.RLN_ORIENT_ROT, randint(-180,180))
   #             transform = rowToAlignment(alignRow, em.ALIGN_2D)
  #              part.setTransform(transform)
 #               trans = part.getTransform()
#                print trans