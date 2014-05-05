# **************************************************************************
# *
# * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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

import os
from os.path import basename

from pyworkflow.em.wizard import EmWizard
from django.shortcuts import render_to_response
from django.http import HttpResponse
from pyworkflow.web.app.em_wizard import *



def wiz_spider_particle_mask_radius(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputParticles
    protParams['label']= "radius"
    protParams['value']= protocol.radius.get()
#    return em_wizard.wiz_particle_mask_radius(protocol, protParams, request)
    return wiz_particle_mask_radius(protocol, protParams, request)


def wiz_spider_particle_mask_radii(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputParticles
    protParams['label']= ["innerRadius", "outerRadius"]
    protParams['value']= [protocol.innerRadius.get(), protocol.outerRadius.get()]
#    return em_wizard.wiz_particles_mask_radii(protocol, protParams, request)
    return wiz_particles_mask_radii(protocol, protParams, request)

def wiz_filter_spider(protocol, request):
    particles = protocol.inputParticles.get()
    
    res = validateParticles(particles) 
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate")
        else:
            context = {'objects': parts,
                       'filterType': protocol.filterType.get(),
                       'filterMode': protocol.filterMode.get(),
                       'usePadding': protocol.usePadding.get(),
                       'protocol': protocol,
                       }
            
            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_filter_spider.html', context)
        
#===============================================================================
# UTILS
#===============================================================================

def get_image_filter_spider(request):
    """
    Function to get the computing image with a spider filter applied
    """
    pars={}
    
    imagePath = request.GET.get('image', None)
    dim = request.GET.get('dim', None)
    filterType = int(request.GET.get('filterType', None))
    pars["filterType"] = filterType
    pars["filterMode"] = int(request.GET.get('filterMode', None))
    pars["usePadding"] = request.GET.get('usePadding', None)    
    pars["op"]="FQ"
    
    # Copy image to filter to Tmp project folder
    outputName = os.path.join("Tmp", "filtered_particle")
    outputPath = outputName + ".spi"
    
    outputLoc = (1, outputPath)
    ih = ImageHandler()
    ih.convert(xmippToLocation(imagePath), outputLoc)
    outputLocSpiStr = locationToSpider(1, outputName)
    
    # check values
    if filterType < 2:
        pars["filterRadius"] = request.GET.get('radius', None)
    else: 
        pars["highFreq"] = float(request.GET.get('highFreq', None))
        pars["lowFreq"] = float(request.GET.get('lowFreq', None))
        if filterType == 2:
            pars["temperature"] = request.GET.get('temperature', None)    

    filter_spider(outputLocSpiStr, outputLocSpiStr, **pars)
    
    # Get output image and update filtered image
    img = xmipp.Image()
    locXmippStr = locationToXmipp(1, outputPath)
    img.read(locXmippStr)
    
    # from PIL import Image
    img = getPILImage(img, dim)
    
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response
