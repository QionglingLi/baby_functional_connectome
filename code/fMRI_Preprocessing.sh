#!/bin/bash
# preprocessing for bcp fmri data
# 12/21/2021 by qiongling
# this script was completed with reference to the MICA pipeline (https://github.com/MICA-MNI/micapipe).

subage=$1
echo $subage

subject="${subage%_*}"
age="${subage#*_}"

Directory=/HeLabData2/qlli/BCPdata/BCP_func_03yr
subDirectory="$Directory"/"$subject"/"$subage"
funcDirectory="$subDirectory"/func
structDirectory="$subDirectory"/anat
processingDirectory="$subDirectory"/processingdirectory
volumetricOutputDirectory="$subDirectory"/volumetricoutputdirectory
templateDirectory=/HeLabData/qlli/atlas/UNC-BCP_4D_Infant_Brain_Volumetric_Atlas_v1
warpDirectory="$subDirectory"/warp
fmap="$subDirectory"/fmap
funcTemDirectory="$funcDirectory"/funtemplate

mkdir $volumetricOutputDirectory $fmap $processingDirectory $funcTemDirectory

# define Modality & template Parameters
warpedstruct=`find "$structDirectory"    -maxdepth 1 -name "*BCP*[[:digit:]]Warped.nii.gz"`
if [[ $warpedstruct == "" ]]; then
	echo "no warp, please do warp first"
	exit
fi

if [ "$age" -le 0 ]; then
Tage="00M"; tage="0Month"
elif [ "$age" -le 1 ]; then
Tage="01M"; tage="1Month"
elif [ "$age" -le 2 ]; then
Tage="02M"; tage="2Month"
elif [ "$age" -le 3 ]; then
Tage="03M"; tage="3Month"
elif [ "$age" -le 4 ]; then
Tage="04M"; tage="4Month"
elif [ "$age" -le 5 ]; then
Tage="05M"; tage="5Month"
elif [ "$age" -le 6 ]; then
Tage="06M"; tage="6Month"
elif [ "$age" -le 7 ]; then
Tage="07M"; tage="7Month"
elif [ "$age" -le 8 ]; then
Tage="08M"; tage="8Month"
elif [ "$age" -le 9 ]; then
Tage="09M"; tage="9Month"
elif [ "$age" -le 10 ]; then
Tage="10M"; tage="10Month"
elif [ "$age" -le 11 ]; then
Tage="11M"; tage="11Month"
elif [ "$age" -le 12 ]; then
Tage="12M"; tage="12Month"
elif [ "$age" -le 15 ]; then
Tage="15M"; tage="15Month"
elif [ "$age" -le 18 ]; then
Tage="18M"; tage="18Month"
elif [ "$age" -le 21 ]; then
Tage="21M"; tage="21Month"
else
Tage="24M"; tage="24Month"
fi

modality="${warpedstruct:0-15:2}"
Modality="$modality""w"

echo "Modality="$Modality
echo "Tage="$Tage

# define templates
hirestemplate=BCP-"$Tage"-"$modality"
lowrestemplate=BCP-"$Tage"-"$modality"_2mm
echo "high resolution template="$hirestemplate
echo "low resolution template="$lowrestemplate

# input func & anatomical MRI
sMRI="$structDirectory"/"$subage"_"$Modality"_nucorrected.nii.gz
fslmaths "$structDirectory"/"$subage"_"$Modality"_brain.nii.gz -bin "$structDirectory"/"$subage"_"$Modality"_mask.nii.gz
fMRI_AP=`find "$funcDirectory"        -maxdepth 1 -name "*rfMRI_REST_AP_[[:digit:]]*.nii.gz"`
SBRef_AP=`find "$funcDirectory"       -maxdepth 1 -name "*rfMRI_REST_AP_SBRef*.nii.gz"`
FieldMap_AP=`find "$funcDirectory"    -maxdepth 1 -name "*SpinEchoFieldMap_AP*.nii.gz"`
FieldMap_PA=`find "$funcDirectory"    -maxdepth 1 -name "*SpinEchoFieldMap_PA*.nii.gz"`

if [[ $fMRI_AP == "" ]]; then
	echo "no fMRI_AP data"
	exit
fi

if [[ $SBRef_AP == "" ]]; then
  echo "no single band image, use first fmri volume"
  fslroi $fMRI_AP "$funcDirectory"/"$subage"_fMRI_AP_1.nii.gz 1 1
  SBRef_AP="$funcDirectory"/"$subage"_fMRI_AP_1.nii.gz 	
 fi

if [[ ! -f "$funcDirectory"/"$subage"_rfMRI_REST_AP_SBRef_reorient.nii.gz ]]; then  
  fslreorient2std $SBRef_AP "$funcDirectory"/"$subage"_rfMRI_REST_AP_SBRef_reorient.nii.gz
fi
SBRef="$funcDirectory"/"$subage"_rfMRI_REST_AP_SBRef_reorient.nii.gz

if [[ ! -f "$funcDirectory"/"$subage"_fMRI_AP_dc.nii.gz ]]; then
	tags=(fMRI_AP FieldMap_AP FieldMap_PA)
	step=-1
	ToProcess="$FieldMap_PA $fMRI_AP $FieldMap_AP"
	for rawdata in ${ToProcess}; do

	  tag=${tags[$step]}
	  step=$(echo $step + 1 | bc)

	  echo "Tag="$tag
	  echo "data="$rawdata


	  # 1. drop first ten TRs and reorient.
	  if [ "$tag" == "fMRI_AP" ]; then
		nifti_tool -cbl -prefix "$funcDirectory"/"$subage"_"$tag"_drop.nii.gz -infiles "$rawdata"'[10..$]'
		fslreorient2std "$funcDirectory"/"$subage"_"$tag"_drop.nii.gz "$funcDirectory"/"$subage"_"$tag"_reorient.nii.gz
	  else
		fslreorient2std "$rawdata" "$funcDirectory"/"$subage"_"$tag"_reorient.nii.gz
	  fi


	  # 2. motion correction with SBRef
	  fslmaths "$funcDirectory"/"$subage"_"$tag"_reorient.nii.gz -Tmean "$funcDirectory"/"$subage"_"$tag"_reorientMean.nii.gz
	  3dvolreg -Fourier -twopass -base $SBRef -zpad 4 -prefix "$funcDirectory"/"$subage"_"$tag"_mc.nii.gz -1Dfile "$volumetricOutputDirectory"/"$subage"_"$tag".1D "$funcDirectory"/"$subage"_"$tag"_reorient.nii.gz
	  fslmaths "$funcDirectory"/"$subage"_"$tag"_mc.nii.gz -Tmean "$funcDirectory"/"$subage"_"$tag"_mcMean.nii.gz
	done
	mv "$funcDirectory"/"$subage"_FieldMap*.nii.gz "$fmap"
  # update SBRef with meanfmri (only for those corrected to first fmri)
	SBRef_AP="$funcDirectory"/"$subage"_"$tag"_mcMean.nii.gz

	fsl_motion_outliers -i "$funcDirectory"/"$subage"_fMRI_AP_reorient.nii.gz -o "$volumetricOutputDirectory"/"$subage"_fMRI_AP_spikeRegressors_FD.1D -s "$volumetricOutputDirectory"/"$subage"_fMRI_AP_metric_FD.1D --fd --thresh=0.5


	# 3. Only do distortion correction if fieldmaps were provided, if not then rename the scan to distortionCorrected
	if [[ $FieldMap_AP == "" ]] || [[ $FieldMap_PA == "" ]]; then
	  echo "shipped topup!"
	  cp $SBRef "$funcDirectory"/"$subage"_SBRef_AP_dc.nii.gz
	  cp "$funcDirectory"/"$subage"_fMRI_AP_mc.nii.gz "$funcDirectory"/"$subage"_fMRI_AP_dc.nii.gz
	else
	  FieldMap_AP=`find "$fmap"        	-maxdepth 1 -name "*FieldMap_AP*_mc.nii.gz"`
	  FieldMap_PA=`find "$fmap"     		-maxdepth 1 -name "*FieldMap_PA*_mc.nii.gz"`
	  FieldMap_AP_mean=`find "$fmap"     	-maxdepth 1 -name "*FieldMap_AP*_mcMean.nii.gz"`
	  FieldMap_PA_mean=`find "$fmap"     	-maxdepth 1 -name "*FieldMap_PA*_mcMean.nii.gz"`

	  flirt -in "$FieldMap_AP_mean" -ref $SBRef -omat "$fmap"/fmap2SBRef_AP.omat
	  flirt -in "$FieldMap_PA_mean" -ref $SBRef -omat "$fmap"/fmap2SBRef_PA.omat
		
	  flirt -in "$FieldMap_AP" -ref $SBRef -applyxfm -init "$fmap"/fmap2SBRef_AP.omat -out "$fmap"/FieldMap_AP_aligned.nii.gz
	  flirt -in "$FieldMap_PA" -ref $SBRef -applyxfm -init "$fmap"/fmap2SBRef_PA.omat -out "$fmap"/FieldMap_PA_aligned.nii.gz
		
	  fslmaths "$fmap"/FieldMap_AP_aligned.nii.gz -Tmean "$fmap"/FieldMap_AP_alignedMean.nii.gz
	  fslmaths "$fmap"/FieldMap_PA_aligned.nii.gz -Tmean "$fmap"/FieldMap_PA_alignedMean.nii.gz
		
	  # distortion correction
	  fslmerge -t "$fmap"/mergeForTopUp.nii.gz "$fmap"/FieldMap_AP_alignedMean.nii.gz "$fmap"/FieldMap_PA_alignedMean.nii.gz
	  topup --imain="$fmap"/mergeForTopUp.nii.gz --datain="$Directory"/topupDataIn.txt --config=b02b0.cnf --out="$fmap"/fmri_topup
	  applytopup --imain="$SBRef" --inindex=1 --datain="$Directory"/topupDataIn.txt --topup="$fmap"/fmri_topup --method=jac --out="$funcDirectory"/"$subage"_SBRef_AP_dc.nii.gz
	  applytopup --imain="$funcDirectory"/"$subage"_fMRI_AP_mc.nii.gz --inindex=1 --datain="$Directory"/topupDataIn.txt --topup="$fmap"/fmri_topup --method=jac --out="$funcDirectory"/"$subage"_fMRI_AP_dc.nii.gz
	fi
fi

hislice=`PrintHeader "$funcDirectory"/"$subage"_fMRI_AP_drop.nii.gz | grep Dimens | cut -d ',' -f 4 | cut -d ']' -f 1 | cut -d ' ' -f 2`
tr=`PrintHeader "$funcDirectory"/"$subage"_fMRI_AP_drop.nii.gz | grep "Voxel Spac" | cut -d ',' -f 4 | cut -d ']' -f 1 | cut -d ' ' -f 2`


# 4. EPI to anatomical registration: distortion corrected SBRef ("$subage"_SBRef_AP_dc) to T2/1
if [[ ! -f "$warpDirectory"/SBRef2"$Modality"_AP.omat ]]; then
flirt -in "$funcDirectory"/"$subage"_SBRef_AP_dc.nii.gz -ref "$sMRI" -omat "$warpDirectory"/SBRef2"$Modality"_AP.omat -dof 6 -out "$funcDirectory"/"$subage"_SBRef_AP_structspace.nii.gz
#flirt -in "$funcDirectory"/"$subage"_fMRI_AP_dc.nii.gz -ref "$sMRI" -applyxfm -init "$warpDirectory"/SBRef2"$Modality"_AP.omat -out "$funcDirectory"/"$subage"_fMRI_AP_structspace.nii.gz
fi


# 5. EPI to BCP age-specific template
if [[ ! -f "$warpDirectory"/"$subage"_BCP-"$Tage"-"$modality"_CollapsedWarp.nii.gz ]]; then
	# combine the ants' displacement field and affine matrix into a single concatenated transformation stored as a displacement field
	echo "combine the ants' displacement field and affine matrix into a single concatenated transformation stored as a displacement field"
	antsApplyTransforms -d 3 -o ["$warpDirectory"/"$subage"_BCP-"$Tage"-"$modality"_CollapsedWarp.nii.gz,1] \
						-t "$warpDirectory"/"$subage"_BCP-"$Tage"-"$modality"1Warp.nii.gz \
						-t "$warpDirectory"/"$subage"_BCP-"$Tage"-"$modality"0GenericAffine.mat \
						-r "$templateDirectory"/"$tage"/"$hirestemplate".nii.gz
fi

if [[ ! -f "$funcDirectory"/"$subage"_fMRI_AP_"$lowrestemplate".nii.gz ]]; then
	# split fmri
	rm -rf "$processingDirectory"/splitfmri "$processingDirectory"/splitfmriwarped "$processingDirectory"/splitbrain "$processingDirectory"/splitstruct
	splitfmri="$processingDirectory"/splitfmri
	splitstruct="$processingDirectory"/splitstruct
	splitbrain="$processingDirectory"/splitbrain
	splitfmriwarped="$processingDirectory"/splitfmriwarped
	mkdir "$splitfmri" "$splitstruct" "$splitfmriwarped" "$splitbrain"

	echo "spliting""$subage"_fMRI_AP_dc.nii.gz
	fslsplit "$funcDirectory"/"$subage"_fMRI_AP_dc.nii.gz "$splitfmri"/"$subage"_fMRI_AP_dc -t

	ind=0000
	while [ "$ind" -lt "$hislice" ]; do
	echo "doing linear registration to native structure..."
	flirt -in "$splitfmri"/"$subage"_fMRI_AP_dc"$ind".nii.gz -ref "$sMRI" -applyxfm -init "$warpDirectory"/SBRef2"$Modality"_AP.omat -out "$splitstruct"/"$subage"_fMRI_AP_structspace"$ind".nii.gz
 
	# extract brain
	fslmaths "$splitstruct"/"$subage"_fMRI_AP_structspace"$ind".nii.gz -mul "$structDirectory"/"$subage"_"$Modality"_mask.nii.gz "$splitbrain"/"$subage"_fMRI_AP_structspace_brain"$ind".nii.gz

	# apply combined displacement field to func
	echo "applying warp to template space"
	antsApplyTransforms -d 3 -n BSpline -v 1 \
						-i "$splitbrain"/"$subage"_fMRI_AP_structspace_brain"$ind".nii.gz \
						-r "$templateDirectory"/"$tage"/"$lowrestemplate".nii.gz \
						-t "$warpDirectory"/"$subage"_BCP-"$Tage"-"$modality"_CollapsedWarp.nii.gz \
						-o "$splitfmriwarped"/"$subage"_fMRI_AP_"$lowrestemplate""$ind".nii.gz
	ind=$(echo $ind + 1 | bc)					
	ind=`echo $ind | awk '{printf("%04d",$0)}'`
	done

	echo "merging""$subage""_fMRI_AP_structspace_brain.nii.gz"
	fslmerge -tr "$splitbrain"/"$subage"_fMRI_AP_structspace_brain.nii.gz "$splitbrain"/"$subage"_fMRI_AP_structspace_brain*.nii.gz "$tr"
	fslmerge -tr "$splitfmriwarped"/"$subage"_fMRI_AP_"$lowrestemplate".nii.gz "$splitfmriwarped"/"$subage"_fMRI_AP_"$lowrestemplate"*.nii.gz "$tr"
	mv "$splitbrain"/"$subage"_fMRI_AP_structspace_brain.nii.gz "$funcDirectory"/"$subage"_fMRI_AP_structspace_brain.nii.gz
	mv "$splitfmriwarped"/"$subage"_fMRI_AP_"$lowrestemplate".nii.gz "$funcDirectory"/"$subage"_fMRI_AP_"$lowrestemplate".nii.gz
fi



# 5_1. registration to common template (BCP-6month)
if [[ ! -f "$funcTemDirectory"/"$subage"_fMRI_AP_template.nii.gz ]]; then
	echo "doing linear registration to common template..."
	flirt -in "$funcDirectory"/"$subage"_fMRI_AP_"$lowrestemplate".nii.gz -ref "$templateDirectory"/6Month/BCP-06M-T1.nii.gz -applyxfm -init "$templateDirectory"/"$tage"/BCP_"$tage"to6Month.mat -out "$funcTemDirectory"/"$subage"_fMRI_AP_template.nii.gz -applyisoxfm 2.0
else
	echo "---------------------------------------------------"
	echo "!!!!  SKIP REGISTRATION TO COMMOM TEMPLATE  !!!!"
	echo "---------------------------------------------------"
fi

# 6_1. smoothing
# spatial smoothing
if [[ ! -f "$funcTemDirectory"/"$subage"_fMRI_AP_smoothed.nii.gz ]]; then
  echo "doing spatial smoothing..."
  3dmerge -1blur_fwhm 4.0 -doall -prefix "$funcTemDirectory"/"$subage"_fMRI_AP_smoothed.nii.gz "$funcTemDirectory"/"$subage"_fMRI_AP_template.nii.gz
else
	echo "---------------------------------------------------"
	echo "!!!!  SKIPPIN SMOOTH  !!!!"
	echo "---------------------------------------------------"
fi

# 7_1. detrend & nuisance regressing
 
# get nuisance regressors (24 motion par + CSF + WM)
mkdir "$volumetricOutputDirectory"/std
if [[ ! -f "$volumetricOutputDirectory"/"$subage"_fMRI_AP_metric_FD.1D ]]; then
  echo "geting outliers..."
  fsl_motion_outliers -i "$funcDirectory"/"$subage"_fMRI_AP_reorient.nii.gz -o "$volumetricOutputDirectory"/"$subage"_fMRI_AP_spikeRegressors_FD.1D -s "$volumetricOutputDirectory"/"$subage"_fMRI_AP_metric_FD.1D --fd --thresh=0.5
fi

if [[ ! -f "$volumetricOutputDirectory"/std/"$subage"_NuisanceRegressors.txt ]]; then
  echo "computing nuiscance regressors..."
  fslmeants -i "$funcTemDirectory"/"$subage"_fMRI_AP_smoothed.nii.gz -o "$volumetricOutputDirectory"/std/"$subage"_WMsignal.mat -m "$templateDirectory"/6Month/BCP-06M-Seg-WM_2mm.nii.gz
  fslmeants -i "$funcTemDirectory"/"$subage"_fMRI_AP_smoothed.nii.gz -o "$volumetricOutputDirectory"/std/"$subage"_CSFsignal.mat -m "$templateDirectory"/6Month/BCP-06M-Seg-CSF_2mm.nii.gz
  fslmeants -i "$funcTemDirectory"/"$subage"_fMRI_AP_smoothed.nii.gz -o "$volumetricOutputDirectory"/std/"$subage"_globalsignal.mat -m "$templateDirectory"/6Month/BCP-06M-BrainMask_2mm.nii.gz
  echo "computing nuiscance regressors by matlab..."
  /usr/local/MATLAB/R2013a/bin/matlab -nodisplay -nodesktop -singleCompThread -r "ComputeRegressors('${subage}','"volumetricoutputdirectory/std"'); quit"
else
	echo "---------------------------------------------------"
	echo "!!!!  SKIP REGRESSORS COMPTUTING   !!!!"
	echo "---------------------------------------------------"
fi

# regress out 24 motion par + CSF + WM
if [[ ! -f "$funcTemDirectory"/"$subage"_fMRI_AP_covregress.nii.gz ]]; then
  echo "doing covariate regressing..."
  fsl_glm -i "$funcTemDirectory"/"$subage"_fMRI_AP_smoothed.nii.gz -d "$volumetricOutputDirectory"/std/"$subage"_NuisanceRegressors.txt --out_res="$funcTemDirectory"/"$subage"_fMRI_AP_covregress.nii.gz
else
	echo "---------------------------------------------------"
	echo "!!!!  SKIP COVARIATE REGRESSING  !!!!"
	echo "---------------------------------------------------"
fi


# 8_1. bandpass filtering
if [[ ! -f "$funcTemDirectory"/"$subage"_fMRI_AP_filtered.nii.gz ]]; then
  echo "doing bandpass filtering..."
  3dTproject -input "$funcTemDirectory"/"$subage"_fMRI_AP_covregress.nii.gz -prefix "$funcTemDirectory"/"$subage"_fMRI_AP_filtered.nii.gz -passband 0.01 0.1
else
	echo "---------------------------------------------------"
	echo "!!!!  SKIP BANDPASS FILTERING  !!!!"
	echo "---------------------------------------------------"
fi

# 9_1. scrubbing
if [[ ! -f "$funcTemDirectory"/"$subage"_fMRI_AP_scrubbed.nii.gz ]]; then
echo "doing scrubbing..." 
  cp "$templateDirectory"/6Month/BCP-06M-BrainMask_2mm.nii.gz "$processingDirectory"/BCP-06M-BrainMask_2mm.nii.gz

  /usr/local/MATLAB/R2013a/bin/matlab -nodisplay -nodesktop -singleCompThread -r "y_Scrubbing('"$funcTemDirectory"/"$subage"_fMRI_AP_filtered.nii.gz', '"$funcTemDirectory"/"$subage"_fMRI_AP_scrubbed.nii.gz', '"$processingDirectory"/BCP-06M-BrainMask_2mm.nii.gz', '"$volumetricOutputDirectory"/"$subage"_TemporalMask.txt', 'linear', 'AfterFiltering', 0, '', '"$tr"'); quit"

  gzip "$funcTemDirectory"/*scrubbed*.nii
  rm -f "$processingDirectory"/BCP-*BrainMask*.nii.gz
else
	echo "---------------------------------------------------"
	echo "!!!!  SKIP SCRUBBING  !!!!"
	echo "---------------------------------------------------"
fi
