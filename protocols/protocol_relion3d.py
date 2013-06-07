#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

from os.path import join, exists
from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
from xmipp import MetaData, Image, MDL_IMAGE, MDL_ITER, MDL_LL, AGGR_SUM, MDL_REF3D, MDL_WEIGHT, \
getBlocksInMetaDataFile, MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDValueEQ, MDL_SAMPLINGRATE
from protlib_utils import runShowJ, getListFromVector, getListFromRangeString
from protlib_parser import ProtocolParser
from protlib_xmipp import redStr, cyanStr
from protlib_gui_ext import showWarning
from protlib_filesystem import xmippExists, findAcquisitionInfo, moveFile
from protocol_ml2d import lastIteration


"""
std::string RelionMainWindow::getCommandLine()
{
    std::string cline;
    if (nr_mpi.getValue() > 1)
        cline="`which relion_refine_mpi`";
    else
        cline="`which relion_refine`";

    // I/O
    // Save the real output name (could be with _ctX for continuation)
    // This name will also be used for the stderr and stdout outputs and the submit script and gui settings filenames
    cline += " --o " + output_name;
    if (is_continue)
    {
        cline += " --continue " + fn_cont.getValue();
    }
    else
    {
        cline += " --i " + fn_img.getValue();
        cline += " --particle_diameter " + floatToString(particle_diameter.getValue());
        cline += " --angpix " + floatToString(angpix.getValue());
        if (run_type == RUN3D)
        {
            cline += " --ref " + fn_ref.getValue();
            if (textToFloat(padding.getValue().c_str()) < 3)
            {
                cline += " --lin --pad " + padding.getValue();
            }
            if (!ref_correct_greyscale.getValue())
            {
                cline += " --firstiter_cc";
            }
            if (ini_high.getValue() > 0.)
            {
                cline += " --ini_high " + floatToString(ini_high.getValue());
            }
        }
    }

    // Optimisation
    cline += " --iter " + floatToString(nr_iter.getValue());
    cline += " --tau2_fudge " + floatToString(tau_fudge.getValue());

    if (do_solvent_flattening.getValue())
    {
        cline += " --flatten_solvent";
        if (fn_mask.getValue().length() > 0)
            cline += " --solvent_mask " + fn_mask.getValue();
    }

    if (do_scale_correction.getValue())
    {
        cline += " --scale";
    }

    if (!is_continue)
    {
        // Only split data in two halves from the start
        if (do_split_random_halves.getValue())
        {
            cline += " --split_random_halves";
        }
    }
    else
    {
        // Only join in a continue run
        if (do_join_random_halves.getValue())
        {
            cline += " --join_random_halves";
        }
    }

    if (!is_continue)
    {
        if (do_ctf_correction.getValue())
        {
            cline += " --ctf";
            if (run_type == RUN3D)
            {
                if (ctf_corrected_ref.getValue())
                    cline += " --ctf_corrected_ref";
            }
            if (only_flip_phases.getValue())
                cline += " --only_flip_phases";
            if (ctf_phase_flipped.getValue())
                cline += " --ctf_phase_flipped";
            if (ctf_intact_first_peak.getValue())
                cline += " --ctf_intact_first_peak";
        }
    }

    // Sampling stuff
    if (!is_continue)
    {
        if (run_type == RUN3D)
        {
            if (strcmp((sym_group.getValue()).c_str(), "O") == 0 ||
                strcmp((sym_group.getValue()).c_str(), "T") == 0 )
                cline += " --sym " + sym_group.getValue();
            else
                cline += " --sym " + sym_group.getValue() + floatToString(sym_nr.getValue());
        }
        cline += " --K " + floatToString(nr_classes.getValue());
    }
    int iover = (int)textToFloat((oversampling.getValue()).c_str());
    if (run_type == RUN3D)
    {
        for (int i = 0; i < 10; i++)
        {
            if (strcmp((sampling.getValue()).c_str(), sampling_options[i].label()) == 0)
            {
                // The sampling given in the GUI will be the oversampled one!
                cline += " --healpix_order " + floatToString((float)i + 1 - iover);
                break;
            }
        }

        if (do_local_ang_searches.getValue())
            cline += " --sigma_ang " + floatToString(sigma_angles.getValue() / 3.);

    }
    if (run_type == RUN2D)
    {
        // The sampling given in the GUI will be the oversampled one!
        cline += " --psi_step " + floatToString(psi_sampling.getValue() * pow(2., iover));
    }
    cline += " --oversampling " + oversampling.getValue();
    cline += " --offset_range " + floatToString(offset_range.getValue());
    // The sampling given in the GUI will be the oversampled one!
    cline += " --offset_step " + floatToString(offset_step.getValue() * pow(2., iover));

    // Always do norm correction
    cline += " --norm";

    cline += " " + other_args.getValue();

    // Running stuff
    cline += " --j " + floatToString(nr_threads.getValue());

    return cline;

}
"""

class ProtRelion3D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.relion3d.name, scriptname, project)
        self.Import = 'from protocol_ml3d import *'
        self.ParamsStr = ''
        if self.DoMlf:
            self.progId = 'mlf'
        else:
            self.progId = 'ml'
        self.ParamsDict['ProgId'] = self.progId
        
        self.addParam('ORoot', self.WorkingDir + '/')
        self.addParam('Extra2D', self.workingDirPath("%s2dextra" % self.progId))
        acquisionInfo = self.findAcquisitionInfo(self.ImgMd)

        if not acquisionInfo is None: 
            md = MetaData(acquisionInfo)
            self.addParam('SamplingRate', md.getValue(MDL_SAMPLINGRATE, md.firstObject()))
                        
    def createFilenameTemplates(self):
        extraRefs = '%(Extra2D)s/iter%(iter)03d/iter_classes.xmd'
        return {
                'iter_logs': "info@" + extraRefs,
                'iter_refs': "classes@" + extraRefs,
                'iter_vols': '%(ExtraDir)s/iter%(iter)03d/iter_volumes.xmd',
                #'refs': 'iter%(iter)06d@' + extraRefs,
                'initial_vols':  '%(ExtraDir)s/initial_volumes.stk',
                'corrected_vols': '%(ExtraDir)s/corrected_volumes.stk',
                'filtered_vols': '%(ExtraDir)s/filtered_volumes.stk',
                'generated_vols': '%(ExtraDir)s/generated_volumes.stk'
                }
        
    def summary(self):
        md = MetaData(self.ImgMd)
        lines = ["Input images:  [%s] (<%u>)" % (self.ImgMd, md.size())]
       
        return lines
    
    def validate(self):
        errors = []
        # TODO: Check if relion is installed
        md = MetaData(self.ImgMd)
        if md.containsLabel(MDL_IMAGE):
            # Check that have same size as images:
            from protlib_xmipp import validateInputSize
            mdRef = MetaData(self.RefMd)
            references = mdRef.getColumnValues(MDL_IMAGE)
            for ref in references:
                if not xmippExists(ref):
                    errors.append("Reference: <%s> doesn't exists" % ref)
            if len(errors) == 0:
                validateInputSize(references, self.ImgMd, errors)
        else:
            errors.append("Input metadata <%s> doesn't contains image label" % self.ImgMd)
            
        return errors 
    
    def defineSteps(self):        
        restart = False
        if restart:            #Not yet implemented
            pass
#        # Restarting a previous run...
#        else:
#            # Execute protocol in the working directory
#            os.chdir(self.WorkingDir)
#            self.restart_MLrefine3D(RestartIter)
        else:
            initVols = self.ParamsDict['InitialVols'] = self.getFilename('initial_vols')
            self.mdVols = MetaData(self.RefMd)
            
            
            self.insertStep('createDir', [self.ExtraDir], path=self.ExtraDir)
            
            self.insertStep('copyVolumes', [initVols], 
                               inputMd=self.RefMd, outputStack=initVols)
            
            if self.DoCorrectGreyScale:
                self.insertCorrectGreyScaleSteps()
                
            if self.DoLowPassFilterReference:
                self.insertFilterStep()
                
            if self.NumberOfReferences > 1:
                self.insertGenerateRefSteps()
                
            
            self.insertML3DStep(self.ImgMd, self.ORoot, self.ParamsDict['InitialVols'], 
                                    self.NumberOfIterations, self.SeedsAreAmplitudeCorrected)
            
            self.insertStep('renameOutput', WorkingDir=self.WorkingDir, ProgId=self.progId)
            
    # Insert the step of launch some program
    def insertRunJob(self, prog, vf=[], files=None, useProcs=True):
        if files is None:
            files = [self.ParamsDict[k] for k in vf]
        NumberOfMpi = NumberOfThreads = 1
        if useProcs:
            NumberOfMpi = self.NumberOfMpi
            NumberOfThreads = self.NumberOfThreads
        self.insertStep('runJob', programname=prog, 
                        params = self.ParamsStr % self.ParamsDict, 
                        verifyfiles = files, 
                        NumberOfMpi = NumberOfMpi, 
                        NumberOfThreads = NumberOfThreads)
    
    def insertML3DStep(self, inputImg, oRoot, initialVols, numberOfIters, amplitudCorrected):
        self.ParamsDict.update({
                         '_ImgMd': inputImg,
                         '_ORoot': oRoot,
                         '_InitialVols': initialVols,
                         '_NumberOfIterations': numberOfIters       
                                })
        self.ParamsStr = "-i %(_ImgMd)s --oroot %(_ORoot)s --ref %(_InitialVols)s --iter %(_NumberOfIterations)d " + \
                         "--sym %(Symmetry)s --ang %(AngularSampling)s %(ExtraParams)s"
#        if self.NumberOfReferences > 1:
#            self.ParamsStr += " --nref %(NumberOfReferences)s"
        if self.NumberOfThreads > 1:
            self.ParamsStr += " --thr %(NumberOfThreads)d"
        if self.DoNorm:
            self.ParamsStr += " --norm"
        
        if self.DoMlf:
            if not self.DoCorrectAmplitudes:
                self.ParamsStr += " --no_ctf"
            if not self.ImagesArePhaseFlipped:
                self.ParamsStr += " --not_phase_flipped"
            if not amplitudCorrected:
                self.ParamsStr += " --ctf_affected_refs"
            if self.HighResLimit > 0:
                self.ParamsStr += " --limit_resolution 0 %(HighResLimit)f"
            self.ParamsStr += ' --sampling_rate %(SamplingRate)f'

        self.ParamsStr += " --recons %(ReconstructionMethod)s "
        
        if self.ReconstructionMethod == 'wslART':
            self.ParamsStr += " %(ARTExtraParams)s"
        else:
            self.ParamsStr += " %(FourierExtraParams)s" 
            
        self.insertRunJob('xmipp_%s_refine3d' % self.progId, [])
        
    def setVisualizeIterations(self):
        '''Validate and set the set of iterations to visualize.
        If not set is selected, only use the last iteration'''
        self.lastIter = lastIter = lastIteration(self)
        self.VisualizeIter = self.parser.getTkValue('VisualizeIter')
        
        if self.VisualizeIter == 'last':
            self.visualizeIters = [lastIter]
        elif self.VisualizeIter == 'all':
            self.visualizeIters = range(1, lastIter+1)
        elif self.VisualizeIter == 'selection':
            selection = getListFromRangeString(self.parser.getTkValue('SelectedIters'))
            self.visualizeIters = [it for it in selection if (it > 0 and it <= lastIter)]
            invalidIters = [it for it in selection if (it <= 0 or it > lastIter)]
            if len(invalidIters):
                print cyanStr("Following iterations are invalid: %s" % str(invalidIters))
        
    def visualize(self):
        self.setVisualizeIterations()
        self.xplotter = None
        for k, v in self.ParamsDict.iteritems():
            if self.parser.hasVariable(k):
                var = self.parser.getVariable(k)
                if var.isVisualize() and v and var.satisfiesCondition():
                    self._visualizeVar(k)
        self.showPlot()
            
    def visualizeVar(self, varName):
        self.xplotter = None
        self.setVisualizeIterations()
        self._visualizeVar(varName)
        self.showPlot()
        
    def _visualizeVar(self, varName):        
        if len(self.visualizeIters) > 0:
            iter = self.visualizeIters[0]
            inputVolVar = {'VisualizeCRVolume': self.getFilename('corrected_vols'), 
                               'VisualizeFRVolume': self.getFilename('filtered_vols'),
                               'VisualizeGSVolume':self.getFilename( 'generated_vols')
                            }
            iterShowVar = { 'VisualizeML3DAvgs': 'iter_refs', 'VisualizeML3DReferences': 'iter_vols'}                          

            if varName in inputVolVar:
                runShowJ(inputVolVar[varName])
            elif varName == 'DoShowStats':
                from protocol_ml2d import launchML2DPlots
                xplotter = launchML2DPlots(self, ['DoShowLL', 'DoShowPmax'])
                self.drawPlot(xplotter)
            elif varName == 'VisualizeClassDistribution':
                for it in self.visualizeIters:
                    self.plotClassDistribution(it)
            elif varName == 'VisualizeAngDistribution':
                for it in self.visualizeIters:
                    self.plotAngularDistribution(it)
            elif varName in iterShowVar:
                for it in self.visualizeIters:
                    runShowJ(self.getFilename(iterShowVar[varName], iter=it))
                    
    def getRefsMd(self, iteration):
        ''' Read the references metadata for a give iteration. '''
        fn = self.getFilename('iter_refs', iter=iteration)        
        return MetaData(fn)
        
    def drawPlot(self, xplotter):
        self.xplotter = xplotter
        
    def showPlot(self):
        if self.xplotter is not None:
            self.xplotter.show()
            
    def plotClassDistribution(self, iteration):
        from numpy import arange
        from protlib_gui_figure import XmippPlotter
        from matplotlib.ticker import FormatStrFormatter
        
        xplotter = XmippPlotter(1, 1, figsize=(4,4),
                                windowTitle="Images distribution - iteration %d" % iteration)
        md = self.getRefsMd(iteration)
        md2 = MetaData()    
        md2.aggregate(md, AGGR_SUM, MDL_REF3D, MDL_WEIGHT, MDL_WEIGHT)
        weights = [md2.getValue(MDL_WEIGHT, objId) for objId in md2]
        nrefs = len(weights)
        refs3d = arange(1, nrefs + 1)
        width = 0.85
        a = xplotter.createSubPlot('3D references weights on last iteration', 'references', 'weight')
        a.set_xticks(refs3d + 0.45)
        a.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
        a.set_xlim([0.8, nrefs + 1])
        a.bar(refs3d, weights, width, color='b')
        self.drawPlot(xplotter)

    def plotAngularDistribution(self, iteration): 
        from protlib_gui_figure import XmippPlotter
        md = self.getRefsMd(iteration)
        md2 = MetaData()
        md2.aggregate(md, AGGR_SUM, MDL_REF3D, MDL_WEIGHT, MDL_WEIGHT)
        nrefs = md2.size()
        figsize = None
        if nrefs == 1:
            gridsize = [1, 1]
            figsize = (4, 4)
        elif nrefs == 2:
            gridsize = [1, 2]
            figsize = (8, 4)
        else:
            gridsize = [(nrefs+1)/2, 2]
            figsize = (8, 12)

        xplotter = XmippPlotter(*gridsize, figsize=figsize, 
                                windowTitle="Angular distribution - iteration %d" % iteration)
        
        for r in range(1, nrefs+1):
            md2.importObjects(md, MDValueEQ(MDL_REF3D, r))  
            plot_title = 'ref %d' % r
            xplotter.plotAngularDistribution(plot_title, md2)
        self.drawPlot(xplotter)
        
                
''' This function will copy input references into a stack in working directory'''
def copyVolumes(log, inputMd, outputStack):
    from protlib_filesystem import deleteFile
    deleteFile(log, outputStack)
    md = MetaData(inputMd)
    img = Image()
    for i, idx in enumerate(md):
        img.read(md.getValue(MDL_IMAGE, idx))
        img.write('%d@%s' % (i + 1, outputStack))
        
def renameOutput(log, WorkingDir, ProgId):
    ''' Remove ml2d prefix from:
        ml2dclasses.stk, ml2dclasses.xmd and ml2dimages.xmd'''
    prefix = '%s2d' % ProgId
    for f in ['%sclasses.stk', '%sclasses.xmd', '%simages.xmd']:
        f = join(WorkingDir, f % prefix)
        nf = f.replace(prefix, '')
        moveFile(log, f, nf)
