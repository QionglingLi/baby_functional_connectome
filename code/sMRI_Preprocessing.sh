#! /bin/bash
# registration for smri data
# 3/30/2022 by qiongling


subage=$1
#subage="106436_00"
echo $subage

subject="${subage%_*}"
age="${subage#*_}"

subDirectory=/HeLabData/qlli/BCPdata/BCP_func_03yr/"$subject"/"$subage"
structDirectory="$subDirectory"/anat
processingDirectory="$subDirectory"/processingdirectory
warpDirectory="$subDirectory"/warp
templateDirectory=/HeLabData/qlli/atlas/UNC-BCP_4D_Infant_Brain_Volumetric_Atlas_v1
MM="0.8mm"

mkdir $processingDirectory
mkdir $warpDirectory

# define Modality & template & Skull Stripping Parameter
if [ "$age" -le 0 ]; then
Tage="00M"; tage="0Month"; Modality=T2w; modality=T2
p=2 # skull strip para
elif [ "$age" -le 1 ]; then
Tage="01M"; tage="1Month"; Modality=T2w; modality=T2
p=2
elif [ "$age" -le 2 ]; then
Tage="02M"; tage="2Month"; Modality=T2w; modality=T2
p=2
elif [ "$age" -le 3 ]; then
Tage="03M"; tage="3Month"; Modality=T2w; modality=T2
p=2
elif [ "$age" -le 4 ]; then
Tage="04M"; tage="4Month"; Modality=T2w; modality=T2
p=2
elif [ "$age" -le 5 ]; then
Tage="05M"; tage="5Month"; Modality=T2w; modality=T2
p=2
elif [ "$age" -le 6 ]; then
Tage="06M"; tage="6Month"; Modality=T2w; modality=T2
p=2
elif [ "$age" -le 7 ]; then
Tage="07M"; tage="7Month"; Modality=T1w; modality=T1
p=1
elif [ "$age" -le 8 ]; then
Tage="08M"; tage="8Month"; Modality=T1w; modality=T1
p=1
elif [ "$age" -le 9 ]; then
Tage="09M"; tage="9Month"; Modality=T1w; modality=T1
p=1
elif [ "$age" -le 10 ]; then
Tage="10M"; tage="10Month"; Modality=T1w; modality=T1
p=1
elif [ "$age" -le 11 ]; then
Tage="11M"; tage="11Month"; Modality=T1w; modality=T1
p=1
elif [ "$age" -le 12 ]; then
Tage="12M"; tage="12Month"; Modality=T1w; modality=T1
p=1
elif [ "$age" -le 15 ]; then
Tage="15M"; tage="15Month"; Modality=T1w; modality=T1
p=1
elif [ "$age" -le 18 ]; then
Tage="18M"; tage="18Month"; Modality=T1w; modality=T1
p=1
elif [ "$age" -le 21 ]; then
Tage="21M"; tage="21Month"; Modality=T1w; modality=T1
p=1
else
Tage="24M"; tage="24Month"; Modality=T1w; modality=T1
p=1
fi

# verify if there is "$Modality" images, if no then try to use another Modality
if [[ "$Modality" = 'T2w' ]] && [[ `find "$structDirectory" -name "*T2w_[[:digit:]]*.nii.gz" | wc -l` -lt 1 ]]; then # no t2w, then try to use t1w
	Modality=T1w; modality=T1
fi

if [[ "$Modality" = 'T1w' ]] && [[ `find "$structDirectory" -name "*T1w_[[:digit:]]*.nii.gz" | wc -l` -lt 1 ]]; then # no t1w, then try to use t2w
	Modality=T2w; modality=T2
fi

# templates
hirestemplate="$templateDirectory"/"$tage"/BCP-"$Tage"-"$modality".nii.gz

# average two images
if [[ ! -f "$structDirectory"/"$subage"_"$Modality".nii.gz ]]; then
	if [[ "$Modality" = 'T2w' ]] && [[ `find "$structDirectory" -name "*T2w_[[:digit:]]*.nii.gz" | wc -l` -gt 1 ]]; then
		fslmaths "$structDirectory"/"$subage"_"$Modality"_1.nii.gz -add "$structDirectory"/"$subage"_"$Modality"_2.nii.gz -mul 0.5 "$structDirectory"/"$subage"_"$Modality".nii.gz
	elif [[ "$Modality" = 'T1w' ]] && [[ `find "$structDirectory" -name "*T1w_[[:digit:]]*.nii.gz" | wc -l` -gt 1 ]]; then
		fslmaths "$structDirectory"/"$subage"_"$Modality"_1.nii.gz -add "$structDirectory"/"$subage"_"$Modality"_2.nii.gz -mul 0.5 "$structDirectory"/"$subage"_"$Modality".nii.gz
	else
		cp "$structDirectory"/"$subage"_"$Modality"_1.nii.gz "$structDirectory"/"$subage"_"$Modality".nii.gz
	fi
fi

# reorient to fsl space
if [[ ! -f "$structDirectory"/"$subage"_"$Modality"_reorient.nii.gz ]]; then
	fslreorient2std "$structDirectory"/"$subage"_"$Modality".nii.gz "$structDirectory"/"$subage"_"$Modality"_reorient
fi

# bias correction using N3
if [[ ! -f "$structDirectory"/"$subage"_"$Modality"_nucorrected.nii.gz ]]; then
	gunzip "$structDirectory"/"$subage"_"$Modality"_reorient.nii.gz
	nii2mnc "$structDirectory"/"$subage"_"$Modality"_reorient.nii "$processingDirectory"/"$subage"_"$Modality"_reorient.mnc
	nu_correct "$processingDirectory"/"$subage"_"$Modality"_reorient.mnc "$processingDirectory"/"$subage"_"$Modality"_nucorrected.mnc
	mnc2nii "$processingDirectory"/"$subage"_"$Modality"_nucorrected.mnc "$processingDirectory"/"$subage"_"$Modality"_nucorrected.nii
	cp "$processingDirectory"/"$subage"_"$Modality"_nucorrected.nii "$structDirectory"/"$subage"_"$Modality"_nucorrected.nii
	gzip "$structDirectory"/*.nii
fi

# extract brain using skullStripping toolkit by Shifeng 2012
if [[ ! -f "$structDirectory"/"$subage"_"$Modality"_brain.nii.gz ]]; then
	fLABEL "$processingDirectory"/"$subage"_"$Modality"_nucorrected.nii -p $p #0 for adult T1 images(default), 1 is for pediatric T1 images, 2 is for neonate T2 images, 3 for unified
	gzip "$processingDirectory"/"$subage"_"$Modality"_nucorrected-strip.nii
	cp "$processingDirectory"/"$subage"_"$Modality"_nucorrected-strip.nii.gz "$structDirectory"/"$subage"_"$Modality"_brain.nii.gz
  fslmaths "$structDirectory"/"$subage"_"$Modality"_brain.nii.gz -bin "$structDirectory"/"$subage"_"$Modality"_mask.nii.gz
fi

# registration to nihpd_asym_template
#if [[ ! -f "$structDirectory"/"$subage"_"$Modality"_nihpd_"$Tage"_"$modality"_"$MM"Warped.nii.gz ]]; then
	echo "registration to "BCP-"$Tage"-"$modality".nii.gz
	antsRegistrationSyN.sh -d 3 -f "$hirestemplate" \
							-m "$structDirectory"/"$subage"_"$Modality"_brain.nii.gz \
							-o "$structDirectory"/"$subage"_BCP-"$Tage"-"$modality"
	mv "$structDirectory"/"$subage"_*.mat "$warpDirectory"/
	mv "$structDirectory"/"$subage"_*Warp.nii.gz "$warpDirectory"/
#fi
