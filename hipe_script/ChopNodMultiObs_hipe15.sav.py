# 
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2013 Herschel Science Ground Segment Consortium
# 
#  HCSS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
# 
#  HCSS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
# 
#  You should have received a copy of the GNU Lesser General
#  Public License along with HCSS.
#  If not, see <http://www.gnu.org/licenses/>.
# 

"""

----------------------------USER'S INSTRUCTIONS--------------------------------
Available via:
 Pipeline -> PACS -> Spectrometer -> Chopped large range scan SED -> Combine obs

Purpose:
  PACS Spectrometer interactive pipeline processing from level 0 to level 2.5 for 
  ChopNod observations. This script calls the interactive pipeline script for 
  several OBSIDs, supposed to belong to the same target.
  Afterwards, it concatenates the slices in a final global spectrum.

DISCLAIMER:
  This script is to be seen as an example of how you can process a few
  obsids together. Don't start here! There are very good reasons why you should
  prefer to process each obsid individually, inspect some intermediate results 
  and tune some parameters to each observation (e.g. spectral flatfielding). 
  You can always combine the end results of individual reductions into one single 
  spectrum via the function 'concatenateSliced'.

WARNING:
  Unlike most of the other interactive pipeline scripts, this one will not run
  "out of the box". Before you use it, you will need to carefully review the 
  input & output variables and tune them to your wishes (obsids.py, script,
  scriptDir, outputDir, ...)
  You are also encouraged to copy the pipeline script you want to run to another
  location on disk and edit the necessary parameters, and run that version here.
  If you use the default pipeline script from the build without any edition, 
  and if you have placed your observations in a non-default location, it will fail
  to find them.

Note: 
  This file is the reference copy, taken from your current installation of HIPE.
  If you need to edit it, we recommend to copy it to a different location first.

Description:
 - The present script includes the loops over obsid and camera. 
   This is not to be regarded as a batch. Its use should be limited to experienced
   users. Conceptually, it is assumed the obsids belong to the same region of 
   the sky, and that the resulting data can be combined into a single and 
   consistent entity (e.g. an SED).
   Hence the present script assumes that a consistent set of two or three
   SED OBSIDs are given on input : B2B + B2A (and/or/not B3A)
 - "obsid", "camera", "oversample" and "upsample" are handled here, 
   but the variables 'verbose', 'useHsa' and 'updateObservationContext' 
   are set in the underlying ipipe script and must still be adjusted in there 
   prior to run this.
 - by default, the final spectrum is saved in a pool whose name starts with 
   the name of the object as found in the meta data of the last observation
   that is processed. This is defined by the variable "objectName" and can
   of course be modified
 - by default, this script will save all data (intermediate and final, in pools
   or in fits files), to the directory specified by 'outputDir'.
   You can change that back to the default behaviour (all pools to "$HOME/.hcss/lstore")
   by simply removing the statement "poolLocation=outputDir" from the calls to
   saveSlicedCopy

Author:
   PACS ICC  

History : 
  2011-04-26 v1 PR
  2013-04-04 v2 PR - include correct3x3 version of final spectra
	 	KE - revision of documentation and comments
"""
import os, time

##### START #####
# Define the 'multiObs' variable, used by the interactive script 
# to branch between interactive and multi-OBSID processing
multiObs = 1

# saveIndividualObsids : 
# 0 - only the global spectrum is saved
# 1 - the slicedCubes, slicedRebinnedCubes & final spectra are saved for every obsid
# in all cases, the final, combined set of spectra is saved (i.e. for all bands) 
saveIndividualObsids = 1


##############################################
###### You need to edit parameters here ######
##############################################
# The pipeline script you will want to run via this multiobs script will be found in
# [where-ever hipe has been put]/scripts/pacs/scripts/ipipe/spec/
# aka: 
# HCSS_DIR = Configuration.getProperty("var.hcss.dir")
# scriptsDir =  HCSS_DIR + "/scripts/pacs/scripts/ipipe/spec/"
#
# In that script there are some parameters that you may need to change, including:
# -verbose
# -useHsa
# -updateObservationContext
# -where your data are - poolLocation and/or poolName must be specfied if your data 
#   are not in your local store, which is usually the directory [HOME]/.hcss/lstore
# -And the parameters of the pipeline tasks that you, via earlier testing, have determined are 
#   best for your observations
# To run the default script, type:
HCSS_DIR = Configuration.getProperty("var.hcss.dir")
#scriptsDir =  HCSS_DIR + "/scripts/pacs/scripts/ipipe/spec/"
#script = scriptsDir + "ChopNodBackgroundNormalizationRange.py"
# But we recommend you instead do something similar to 

scriptsDir =  "/home/bettyjo/yaolun/programs/misc/hipe_script/"
script = scriptsDir + "ChopNodBackgroundNormalizationRange.py"

# Chose the obsids
obsids = {}
obsids[1342269932] = "SEDA" 
obsids[1342269933] = "SEDB"
"""
# You can also read these from a file, which should look like this
# --- start of file ---
obsids = {}
obsids[1342196694] = "SEDA" 
obsids[1342196695] = "SEDB"
# --- end  of  file ---
# And then you run that file (here we call it obsids.py)
#   - if it is in directory you started HIPE from, and called obsids.py:
obsidsDir = Configuration.getProperty("user.dir")
execfile(obsidsDir + "obsids.py")
#   - if it is in another directory
obsidsDir = "/Users/me/"
execfile(obsidsDir + "obsids.py")
"""
# The few pipeline parameters you can set here:s
# Rebinning Parameters
oversample = 8
upsample   = 7
strovup = "ov"+str(oversample)+"_up"+str(upsample)

######################################################
###### You need to edit parameters here as well ######
######################################################


# Example: outputDir = "/home/me/Herschel/PacsSpectro/pipelineOutput/"
# When saveOutput is True, nameBasis will be used as basis for the filenames of all outputs
# outputDir  = str(Configuration.getWorkDir())+"/pacsSpecOut/"
outputDir = '/home/bettyjo/yaolun/L1251B/ov8_up7/'
if (not os.path.exists(outputDir)): os.mkdir(outputDir)

prefix=""
if outputDir: prefix = outputDir+"/"

# For traceability, you can use the buildNumber in the file or directory name
buildNumber = str(Configuration.getProjectInfo().track) + '.' + str(Configuration.getProjectInfo().build)

# Create file for tracking the progress
trackfilename = prefix+"ChopNodSEDMultiObsTracking.txt"
trackfile = open(trackfilename,'w')
trackfile.write("START \n")
trackfile.close()


# Structure holding the final cubes for every pair [obsid,camera]
finalCubeList = []

starttime = time.time()


############################################################
##### Run the pipeline over all the cameras and obsids #####
############################################################
for obsid in obsids.keys():
	for camera in ["blue","red"]:
		# print outs to keep you up to date with progress
		print
		print "START ", obsid, camera
		print
		#target = obsids[obsid].split("_")[0]
		trackfile = open(trackfilename,'a')
		trackfile.write("START " + str(obsid) + " " + camera + "\n")
		trackfile.close()
		obsidstarttime = time.time()
		#
		# Process the observation
		# -----------------------
		execfile(script)
		#
		# Save the results
		nameBasis = "OBSID_"+str(obsid)+"_"+obsids[obsid]+"_"+camera+"_"+buildNumber.replace('.','_')+"_"+strovup
		if saveIndividualObsids:
			saveSlicedCopy(slicedCubes, nameBasis+"_slicedCubes", poolLocation=outputDir)
			saveSlicedCopy(slicedRebinnedCubes, nameBasis+"_RebinnedCubes", poolLocation=outputDir)
			saveSlicedCopy(slicedFinalCubes, nameBasis+"_FinalCubes", poolLocation=outputDir)
		#
		trackfile = open(trackfilename,'a')
		trackfile.write("END   " + str(obsid) + " " + camera + " Duration: "+ str("%7.1f\n" % (time.time() - obsidstarttime)) +"\n")
		trackfile.close()
		#
		# Gather the final spectra in one single structure, later easier to handle, save & plot
		theseFinalCubes = getSlicedCopy(slicedFinalCubes)
		finalCubeList.append(theseFinalCubes)

# Merge all finalCubes into one single slicedProduct
allFinalCubes = concatenateSliced(finalCubeList)

# save the final result in a single pool
objectName = "NoObjectName"
try:
	objectName = obs.meta["object"].value.replace(" ","").replace("+","plus").replace(".","_")
	objectName += "_"+str(obs.meta["odNumber"].value)
except:
	pass
name = objectName + "_ChopNodMultiObs_allFinalCubes"
saveSlicedCopy(allFinalCubes , name, poolLocation=outputDir)

if verbose:
	slicedSummary(allFinalCubes)
	# Display the final, composite spectrum, made of all slices of all OBSIDs
	x,y = 2,2
	pfinal = plotCubes(allFinalCubes,[],x=x,y=y)


#########################################################
#### IF POINT SOURCE : extract central spaxel
####                   apply point source correction
####                   save to FITS
#fullSpec = Spectrum1d()
#segments = Int1d()
#fullSpec9 = Spectrum1d()
#segments9 = Int1d()
#fullSpecCorr3x3 = Spectrum1d()
#segmentsCorr3x3 = Int1d()
#
#### The spectrum "Corrected3x3YES", and the parameters used to produce it must 
#### be used with care, so please carefully review the parameters below (see the PDRG)
#target = slicedFinalCubes.meta["object"].value.replace(" ","_")
#
#smoothing = 'wavelet'
#nLowFreq            = 4
## or
#smoothing = 'filter'
#gaussianFilterWidth = 50
#medianFilterWidth   = 15
#
#for slice in range(len(slicedFinalCubes.refs)):
#    # a. Extract central spectrum, incl. point source correction (c1)
#    # b. Extract SUM(central_3x3 spaxels), incl. point source correction (c9)
#    # c. Scale central spaxel to the level of the central_3x3 (c129 -> See PDRG & print extractCentralSpectrum.__doc__)
#    central1,central9,central129 = extractCentralSpectrum(slicedFinalCubes.get(slice), smoothing=smoothing, width=gaussianFilterWidth, preFilterWidth=medianFilterWidth, nLowFreq=nLowFreq, calTree=calTree, verbose=verbose)
#    #
#    # Save to Fits
#    if saveIndividualObsids:
#	name = "OBSID_"+str(obsid)+"_"+target+"_"+camera+"_centralSpaxel_PointSourceCorrected_Corrected3x3NO_slice_"
#	simpleFitsWriter(product=central1,file = name+str(slice).zfill(2)+".fits")
#	name = "OBSID_"+str(obsid)+"_"+target+"_"+camera+"_central9Spaxels_PointSourceCorrected_slice_"
#	simpleFitsWriter(product=central9,file = name+str(slice).zfill(2)+".fits")
#	name = "OBSID_"+str(obsid)+"_"+target+"_"+camera+"_centralSpaxel_PointSourceCorrected_Corrected3x3YES_slice_"
#	simpleFitsWriter(product=central129,file = name+str(slice).zfill(2)+".fits")
#    #
#    if verbose:
#        centralSpec = central129.spectrum1d
#        try:
#            openVariable("centralSpec", "Spectrum Explorer")
#        except:
#            print "Spectrum Explorer only works within hipe gui, not when calling hipe from the command line"
#    #
#    # Concatenate all bands (slices) to one spectrum1d, keeping them as separate segments
#    spec = central1.spectrum1d
#    fullSpec.concatenate(spec)
#    segments.append(Int1d(spec.flux.length(), slice))
#    spec9 = central9.spectrum1d
#    fullSpec9.concatenate(spec9)
#    segments9.append(Int1d(spec9.flux.length(), slice))
#    specCorr3x3 = central129.spectrum1d
#    fullSpecCorr3x3.concatenate(specCorr3x3)
#    segmentsCorr3x3.append(Int1d(specCorr3x3.flux.length(), slice))
#
## Include the column segment in the final Spectrum1d & save to fits
#fullSpec.setSegment(segments)
#fullSpec9.setSegment(segments9)
#fullSpecCorr3x3.setSegment(segmentsCorr3x3)
#name = objectName + "_ChopNodSEDMultiObs_1d_centralSpaxel_correct3x3_NO_PointSourceCorrected.fits"
#simpleFitsWriter(SimpleSpectrum(fullSpec),prefix+name)
#name = objectName + "_ChopNodSEDMultiObs_1d_central9Spaxels_PointSourceCorrected.fits"
#simpleFitsWriter(SimpleSpectrum(fullSpec9),prefix+name)
#name = objectName + "_ChopNodSEDMultiObs_1d_centralSpaxel_correct3x3_YES_PointSourceCorrected.fits"
#simpleFitsWriter(SimpleSpectrum(fullSpecCorr3x3),prefix+name)
#
#if verbose:
#    # Display in spectrum explorer
#    try:
#	    openVariable("fullSpec", "Spectrum Explorer")
#	    openVariable("fullSpec9", "Spectrum Explorer")
#	    openVariable("fullSpecCorr3x3", "Spectrum Explorer")
#    except:
#        print "Info: openVariable could not open Spectrum Explorer"

trackfile = open(trackfilename,'a')
trackfile.write("END			Total Duration: "+ str("%7.1f\n" % (time.time() - starttime)) +"\n")
trackfile.close()

