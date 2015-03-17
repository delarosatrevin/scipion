#!/usr/bin/env python

from pyworkflow.em.viewer import ChimeraClient, ChimeraProjectionClient
import os, sys
import argparse
import xmipp

def main():
    commonParser = argparse.ArgumentParser(add_help=False, prog='Chimera Client')
    commonParser.add_argument('--input', help='Volume to visualize', required=True)
    commonParser.add_argument('--samplingRate', help='Volume sampling rate')

    commonParser.add_argument('--angDistFile', help='Angular distribution file')
    commonParser.add_argument('--spheresColor', default='red', help='Angular distribution spheres color')
    commonParser.add_argument('--spheresDistance', type=int, help='Angular distribution spheres distance')
    commonParser.add_argument('--spheresMaxRadius', type=int, help='Angular distribution spheres max radius')


    parentParser = argparse.ArgumentParser(add_help=False, prog='Chimera Client')
    subparsers = parentParser.add_subparsers(dest='cmd')
    viewerParser = subparsers.add_parser('viewer', help='Display volume', parents=[commonParser])

    projectorParser = subparsers.add_parser('projector', help='Projector mode displays volume projection.', parents=[commonParser])
    projectorParser.add_argument('--projectionSize', help='Projection window dimensions', type=int)
    projectorParser.add_argument('--paddingFactor', default=1, type=float, help='Projection padding factor')
    projectorParser.add_argument('--maxFreq', default=0.5, type=float, help='Projection spline degree')
    projectorParser.add_argument('--splineDegree', default='BSPLINE2',
                        choices=['NEAREST', 'LINEAR', 'BSPLINE2', 'BSPLINE3', 'BSPLINE4'], help='projection window dimensions')


    splineDegreeDict = {'NEAREST': xmipp.NEAREST, 'LINEAR': xmipp.LINEAR,
                        'BSPLINE2': xmipp.BSPLINE2, 'BSPLINE3': xmipp.BSPLINE3, 'BSPLINE4': xmipp.BSPLINE4}

    args = parentParser.parse_args()
    #print args
    volfile = args.input
    voxelSize= args.samplingRate if hasattr(args, 'samplingRate') else None
    angularDistFile = args.angDistFile if hasattr(args, 'angDistFile') else None
    spheresColor = args.spheresColor if hasattr(args, 'spheresColor') else None
    spheresDistance = args.spheresDistance if hasattr(args, 'spheresDistance') else None
    spheresMaxRadius = args.spheresMaxRadius if hasattr(args, 'spheresMaxRadius') else None

    if args.cmd == 'viewer':
        ChimeraClient(volfile, angularDistFile=angularDistFile, spheresColor=spheresColor, spheresDistance=spheresDistance, spheresMaxRadius=spheresMaxRadius, voxelSize=voxelSize)
    else:
        projectionSize = args.projectionSize if hasattr(args, 'projectionSize') else None
        splineDegree = splineDegreeDict.get(args.splineDegree)
        paddingFactor= args.paddingFactor
        maxFreq = args.maxFreq
        ChimeraProjectionClient(volfile,
                                angularDistFile=angularDistFile,
                                spheresColor=spheresColor,
                                spheresDistance=spheresDistance,
                                spheresMaxRadius=spheresMaxRadius,
                                size=projectionSize,
                                paddingFactor=paddingFactor,
                                maxFreq=maxFreq,
                                splineDegree=splineDegree,
                                voxelSize=voxelSize)

    
def which(file):
    for path in os.environ["PATH"].split(":"):
        if os.path.exists(path + "/" + file):
                return path + "/" + file

    return None
    
if __name__ == '__main__':
    main()
    
    

