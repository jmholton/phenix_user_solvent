#! /bin/tcsh -f
#
# remove bulk solvent from Fobs  - ala "squeeze"   - James Holton  3-30-26
#
# requires Tc_maskify.com map_func.com map_scaleB.com
#
# default to use refmac outputs
set pdbfile = refmacout.pdb
set mtzfile = refmacout.mtz

# seems best to remove hydrogen before surface calc
set exclude_H = 1

# new Fobs with solvent removed
set outfile = refme_minusol.mtz
# new Fobs with solvent difference features as a partial structure
set refmacoutfile = refme_Fpart.mtz

# labels in mtz files
set FP = "auto"
set FreeR_flag = "auto"
set reso = "auto"
set oversample = 10

# Fcalc to subtract initially
set FC_subtract = FC_ALL_LS,PHIC_ALL_LS
# Fcalc to add back to difference map in the end to form new Fobs
set FC_addback = FC,PHIC

# do more than squash, over-shoot solvent correction
set overshoot_scale = -0.5

# use an external map as a guide for grid, limits and axes
set guide_map = ""

# phenix _f_model.mtz input mode (auto-detected below); set 1 to force
set phenix_mode = 0

# exclude Free-R flagged reflections from difference map
set exclude_freeR = 1

# select only statistically significant difference features
set significance = 1

# B factor for enlarging Shannon voxels
set ShannonB = 9.483
set ShannonB = auto
#set ShannonB = 0

# B factor to use when calculating initial difference map
set fftB = 0
set fftreso = auto

# B factor to use when passing sharp probability mask through fft filter
# setting to zero means no filter, very small value means just back-and-forth fft
set filterB = 1e-6
# number of times to pass mask through fft filter
set recycles = 5

# max number of voxels to grow the mask in one pass (can be faster if smaller)
set max_pix_grow = 3

# straighten the Wilson plot of the fofc map
set wilsonify = 1
# high and low fraction of (sin(theta)/lambda)^2 space to fit to straight lines
set wilsonify_hifrac = 0.2
set wilsonify_lofrac = 0.25

# resolution filter applied to the difference density before masking (built into
# the same per-reflection scale the Wilsonification uses):
#   wilson   = straighten the Wilson plot to the low-res slope (default)
#   bandpass = Gaussian band-pass in stol^2, centered on the signal band (SNR peak):
#              H = exp(-(stol^2-bandpass_center)^2 / (2*bandpass_width^2))
#   blur     = monotonic B-factor rolloff: H = exp(-blur_B*stol^2)
#   PLANNED (data-driven, need per-shell noise power Pnn=<sigF^2> from $mtzfile):
#     wiener = Pss/(Pss+Pnn)   (Pss = <|Fdiff|^2>-Pnn)   -- self-tuning band-pass
#     snr    = Pss/Pnn         (matched filter)
set resfilter = wilson
# band-pass centre/width in stol^2 = (sin(theta)/lambda)^2; 0.05 ~ 2.2 A (the 1.9-2.7A SNR peak)
set bandpass_center = 0.05
set bandpass_width  = 0.05
set blur_B = 20

# pre-scale the FOFC map
set fofc_prescale = 2

# select buffer zone around significant difference peaks
set radius = auto
set bevel = 0.2

# suppress difference features in the protein region
set exclude_protein = 1
# probabilities considered definitely there vs unlikely
# protein map will be range-stretched to make this range equal 0-1
set protein_highprob = 0.9
set protein_lowprob = 0.1

# probabilities considered definitely there vs unlikely
#set highprob = 0.87
#set lowprob = 0.145
set highprob = 0.9
set lowprob = 0.1

# allow selecting only one peak at a time
set nearpeak = 0
set nearpeak_radius = 3
set nearpeak_bevel  = 1
set nearpeak_recycles = 1

# use segmentation to establish borders of peak
set floodfill = 0
set floodfill_retries = 100
set floodfill_sigma = 0.5
set floodfill_Bsmooth = 1
set floodfill_maskBsmooth = 1
set floodfill_rethresh = 0.8
set floodfill_maxradius = 5
set floodfill_fftmargin = 2
set floodfill_filterB = 5
set floodfill_recycles = 3
set floodfill_highprob = 0.9
set floodfill_lowprob = 0.1

# default temporary file location
set tempfile = /dev/shm/${USER}/squish_temp_$$_
set debug = 0

set logfile = details.log

set path = ( $path `dirname $0` )

# read the command line to update variables and other settings
foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $Arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if( $assign ) then
      # re-set any existing variables
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" =~ exclude_rfree* ) then
        set exclude_freeR = "$val"
        continue
      endif
      if("$key" == toppeak ) then
        set nearpeak = "$val"
        continue
      endif
      echo "WARNING: did not understand $Arg"
      continue
    endif
    # no equal sign
    if("$Arg" =~ *.pdb ) set pdbfile = "$Arg"
    if("$Arg" =~ *.mtz ) set mtzfile = "$Arg"
    if("$Arg" =~ *.map ) set guide_map = "$Arg"
    
    if("$Key" == "debug") set debug = "1"
end

# establish temp file location
foreach tempfile ( $tempfile /dev/shm/${USER}/squishtemp_$$_ /tmp/${USER}/squishtemp_$$_ ~/squishtemp_$$_ )
    set tempdir = `dirname $tempfile`
    if(! -w "$tempdir") mkdir -p $tempdir
    if(-w "$tempdir") break
end
if(! -w "$tempdir") then
    set BAD = "cannot make temporary files"
    goto exit
endif
set t = "${tempfile}"
set tempdir = `dirname $tempfile`

# round of things that should be integers
foreach variable ( exclude_H exclude_freeR exclude_protein floodfill_retries recycles nearpeak_recycles floodfill_recycles wilsonify debug )
    set before = `eval echo \$$variable`
    set $variable = `echo $before | awk '{print int($1)}'`
    set after = `eval echo \$$variable`
    if( "$before" != "$after" ) echo "WARNING: $variable $before rounded to $after"
end

if(! -e "$pdbfile") then
    set BAD = "need pdb file: $pdbfile"
    goto exit
endif
if(! -e "$mtzfile") then
    set BAD = "need refmac output mtz file: $mtzfile"
    goto exit
endif

foreach dependency ( map_scaleB.com map_func.com Tc_maskify.com addup_maps_runme.com pick.com )
   echo -n "using: "
   which $dependency
   if( $status ) then
       set BAD = "need $dependency in "'$'"path"
       goto exit
   endif
end

# examine the mtz file
echo header | mtzdump hklin $mtzfile |\
   tee ${t}mtzdump.txt |\
   awk '/^ H K L /{for(i=1;i<=NF;++i)label[i]=$i}\
        /^ H H H /{for(i=1;i<=NF;++i)print $i,label[i]}' |\
   cat >! ${t}labels.txt

# ---- phenix _f_model.mtz input support (see phenix2squish.py) ----
# phenix output has FMODEL (total scaled model) but no refmac FC_ALL_LS.
set has_fmodel  = `awk '$2=="FMODEL"' ${t}labels.txt | wc -l`
set has_fcallls = `awk '$2=="FC_ALL_LS"' ${t}labels.txt | wc -l`
if( $has_fmodel && ! $has_fcallls ) then
    set phenix_mode = 1
    set has_freerflag = `awk '$2=="FreeR_flag"' ${t}labels.txt | wc -l`
    set has_rfree     = `awk '$2=="R_FREE_FLAGS"' ${t}labels.txt | wc -l`
    if( $has_rfree && ! $has_freerflag ) then
        # raw _f_model.mtz: convert (flip free flag to CCP4 conv, add Fpart=K_MASK*FMASK)
        echo "phenix _f_model.mtz detected: converting with phenix2squish.py"
        set adapter = `which phenix2squish.py`
        if( "$adapter" == "" ) then
            set BAD = "phenix input needs phenix2squish.py and phenix.python in "'$'"path"
            goto exit
        endif
        phenix.python $adapter $mtzfile ${t}squishin.mtz >> $logfile
        if( $status || ! -e ${t}squishin.mtz ) then
            set BAD = "phenix2squish.py failed (need phenix.python sourced)"
            goto exit
        endif
        set mtzfile = ${t}squishin.mtz
        # rebuild header/labels from the converted mtz
        echo header | mtzdump hklin $mtzfile |\
           tee ${t}mtzdump.txt |\
           awk '/^ H K L /{for(i=1;i<=NF;++i)label[i]=$i}\
                /^ H H H /{for(i=1;i<=NF;++i)print $i,label[i]}' |\
           cat >! ${t}labels.txt
    else
        echo "phenix-style mtz detected (already converted)"
    endif
    # adopt phenix column vocabulary, only overriding the refmac defaults
    if( "$FP" == "auto" ) set FP = FOBS
    set SIGFP = SIGFOBS
    if( "$FC_subtract" == "FC_ALL_LS,PHIC_ALL_LS" ) set FC_subtract = FMODEL,PHIFMODEL
    if( "$FC_addback"  == "FC,PHIC" ) set FC_addback = FMODEL,PHIFMODEL
    echo "phenix mode: FP=$FP FC_subtract=$FC_subtract FC_addback=$FC_addback"
endif

set SGnum = `awk '/Space group =/{print $NF+0}' ${t}mtzdump.txt | tail -n 1`
set symops = `awk -v n=$SGnum '$1==n{print $3;exit}' ${CLIBD}/symop.lib`
set CELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' ${t}mtzdump.txt`
echo $CELL |\
awk 'NF==6{DTR=atan2(1,1)/45; A=cos(DTR*$4); B=cos(DTR*$5); G=cos(DTR*$6); \
 skew = 1 + 2*A*B*G - A*A - B*B - G*G ; if(skew < 0) skew = -skew;\
 printf "%.3f\n", $1*$2*$3*sqrt(skew)}' |\
cat >! ${t}volume
set CELLvolume = `cat ${t}volume`
rm -f ${t}volume

if( "$reso" == "auto" ) then
   set reso = `awk '/Resolution Range/{getline;getline;print $6}' ${t}mtzdump.txt`
   echo "resolution in $mtzfile is $reso"
endif
if( "$FP" == "auto" ) then
   cat  ${t}labels.txt |\
   egrep -v FC |\
   awk '$1=="F"{print $2}' >! ${t}Fs.txt
   set FP = `head -n 1 ${t}Fs.txt`
   if( "$FP" == "" ) then
      set BAD = "no Fobs in $mtzfile"
      goto exit
   endif
   set SIGFP = `grep "Q SIG$FP" ${t}labels.txt | awk '{print $2;exit}'`
   if( "$SIGFP" == "" ) then
      set SIGFP = `egrep "^Q " ${t}labels.txt | awk '{print $2;exit}'`
   endif
   if( "$SIGFP" == "" ) then
      set BAD = "no sigFobs in $mtzfile"
      goto exit
   endif
   echo "selecting Fobs = $FP $SIGFP"
endif
if( "$FreeR_flag" == "auto" ) then
   cat  ${t}labels.txt |\
   awk '$1=="I" && tolower($2)~/free/{print $2}' >! ${t}free.txt
   set FreeR_flag = `head -n 1 ${t}free.txt`
   if( "$FreeR_flag" == "" ) then
      set FreeR_flag = `egrep "^I " ${t}labels.txt | awk '{print $2;exit}'`
   endif
   if( "$FreeR_flag" == "" ) then
      set BAD = "no FreeR_flag in $mtzfile"
      goto exit
   endif
   echo "selecting FreeR_flag = $FreeR_flag"
endif
set hires = `echo $reso $oversample | awk '{print (($1**3)/$2)**0.3333}'`
echo "setting hires: $hires A for ${oversample}x over-sampling of map grid"

foreach FC ( $FC_subtract $FC_addback )
  echo $FC |\
   awk -F "," '{print $1,$2}' |\
   cat - ${t}labels.txt |\
   awk 'NR==1{F=$1;P=$2;next}\
     {++seen[$2]}\
     END{print seen[F],seen[P]}' >! ${t}test.txt
  set test = `cat ${t}test.txt`
  if( "$test" != "1 1" ) then
    set BAD = "no $FC in $mtzfile"
    goto exit
  endif
end

if( -e "$guide_map") then
   echo "using map grid from $guide_map"
   echo "go" | mapdump mapin $guide_map >! ${t}mapdump.txt
  set mapSG     = `awk '/Space-group/{print $NF; exit}' ${t}mapdump.txt`
  set mapCELL   = `awk '/Cell dimensions/{print $(NF-5), $(NF-4), $(NF-3), $(NF-2), $(NF-1), $NF; exit}' ${t}mapdump.txt`
  set mapGRID   = `awk '/Grid sampling/{print "GRID",$(NF-2), $(NF-1), $NF; exit}' ${t}mapdump.txt`
  set mapLIMITS = `awk '/Fast, medium, slow/{o[$(NF-2)]=1;o[$(NF-1)]=2;o[$NF]=3; print b[o["X"]],e[o["X"]], b[o["Y"]],e[o["Y"]], b[o["Z"]],e[o["Z"]]; exit} /Start and stop points/{b[1]=$(NF-5); e[1]=$(NF-4); b[2]=$(NF-3); e[2]=$(NF-2); b[3]=$(NF-1); e[3]=$NF}' ${t}mapdump.txt`
  set mapAXIS   = `awk '/Fast, medium, slow/{print $(NF-2), $(NF-1), $NF}' ${t}mapdump.txt | awk '! /[^ XYZ]/'`
   echo "$mapGRID"
   echo "xyzlim $mapLIMITS"
   echo "axis $mapAXIS"
else
  set mapGRID = "reso $hires"
  set mapLIMITS = "asu"
  set mapAXIS = "X Y Z"
endif

echo -n "settings: "
cat << EOF | awk 'BEGIN{ORS=" "} {print}'
pdbfile=$pdbfile
mtzfile=$mtzfile
outfile=$outfile
refmacoutfile=$refmacoutfile
FP=$FP
FreeR_flag=$FreeR_flag
reso=$reso
oversample=$oversample
FC_subtract=$FC_subtract
FC_addback=$FC_addback
overshoot_scale=$overshoot_scale
guide_map=$guide_map
exclude_freeR=$exclude_freeR
ShannonB=$ShannonB
filterB=$filterB
wilsonify=$wilsonify
recycles=$recycles
radius=$radius
bevel=$bevel
exclude_protein=$exclude_protein
exclude_H=$exclude_H
protein_highprob=$protein_highprob
protein_lowprob=$protein_lowprob
highprob=$highprob
lowprob=$lowprob
nearpeak=$nearpeak
nearpeak_radius=$nearpeak_radius
nearpeak_bevel=$nearpeak_bevel
floodfill=$floodfill
floodfill_sigma=$floodfill_sigma
floodfill_Bsmooth=$floodfill_Bsmooth
floodfill_maskBsmooth=$floodfill_maskBsmooth
floodfill_rethresh=$floodfill_rethresh
floodfill_maxradius=$floodfill_maxradius
floodfill_recycles=$floodfill_recycles
floodfill_highprob=$floodfill_highprob
floodfill_lowprob=$floodfill_lowprob
tempfile=$tempfile
debug=$debug
EOF
echo ""
echo ""
touch $logfile


echo "making $FP - $FC_subtract difference map"
set FC = `echo $FC_subtract | awk -F "," '{print $1}'`
set PHIC = `echo $FC_subtract | awk -F "," '{print $2}'`
if( $exclude_freeR ) then
   echo "excluding $FreeR_flag hkls from fofc.map calculation"
   set FREE = "FREE=$FreeR_flag"
else
   echo "including all hkls in fofc.map calculation"
   set FREE = ""
endif
set resofft = ""
if( "$fftreso" != "auto" ) set resofft = "reso $fftreso"

fft hklin $mtzfile mapout ${t}ffted.map << EOF >> $logfile
labin F1=$FP F2=$FC PHI=$PHIC $FREE
$mapGRID
$resofft
EOF
mapmask mapin ${t}ffted.map mapout fofc0.map << EOF >> $logfile
xyzlim $mapLIMITS
axis $mapAXIS
EOF
echo "multiplying fofc map by $fofc_prescale B=$fftB"
map_scaleB.com fofc0.map scale=$fofc_prescale B=$fftB output=fofc.map >> $logfile

echo | mapdump mapin fofc.map >! ${t}mapdump.txt
set fofc_stats = `awk '/density/{print $NF}' ${t}mapdump.txt`
set fofc_rmsd = `echo $fofc_stats | awk '{print $4}'`
set mapGRID   = `awk '/Grid sampling/{print "GRID",$(NF-2), $(NF-1), $NF; exit}' ${t}mapdump.txt`
set mapCELL   = `awk '/Cell dimensions/{print $(NF-5), $(NF-4), $(NF-3), $(NF-2), $(NF-1), $NF; exit}' ${t}mapdump.txt`
set voxel_size = `echo $mapGRID $mapCELL | awk '{print $5/$2,$6/$3,$7/$4}' | awk '{print $1,$2,$3,($1*$2*$3)^(1./3)}'`
set voxels = `awk '/Number of columns, rows, sections/{print $7*$8*$9}' ${tempfile}mapdump.txt`
set size = `echo $voxels | awk '{print 4*$1}'`
set mapheadersize = `ls -lL fofc.map | awk -v size=$size '{print $5-size}'`


# diagnostic maps below use refmac-only columns (DELFWT/FWT/FC_ALL*); skip for phenix
if( ! $phenix_mode ) then
fft hklin $mtzfile mapout ${t}ffted.map << EOF >> $logfile
labin F1=$FP PHI=$PHIC
$mapGRID
EOF
mapmask mapin ${t}ffted.map mapout Fobs.map << EOF >> $logfile
xyzlim $mapLIMITS
axis $mapAXIS
EOF

fft hklin $mtzfile mapout ${t}ffted.map << EOF >> $logfile
labin F1=DELFWT PHI=PHDELWT
$mapGRID
EOF
mapmask mapin ${t}ffted.map mapout mfodfc.map << EOF >> $logfile
xyzlim $mapLIMITS
axis $mapAXIS
EOF

fft hklin $mtzfile mapout ${t}ffted.map << EOF >> $logfile
labin F1=FWT PHI=PHWT
$mapGRID
EOF
mapmask mapin ${t}ffted.map mapout 2mfodfc.map << EOF >> $logfile
xyzlim $mapLIMITS
axis $mapAXIS
EOF




foreach C ( C C_ALL C_ALL_LS )
fft hklin $mtzfile mapout ${t}ffted.map << EOF >> $logfile
labin F1=F$C PHI=PHI$C
$mapGRID
EOF
mapmask mapin ${t}ffted.map mapout F${C}.map << EOF >> $logfile
xyzlim $mapLIMITS
axis $mapAXIS
EOF
end
endif
# end of phenix-skipped diagnostic maps





if( ! $exclude_protein ) then
    goto skip_protein
endif

if( $exclude_H ) then
    cat $pdbfile |\
    awk '! /^ATOM|^HETAT/{print;next}\
      {Ee = substr($0,77,4);gsub(" ","",Ee)}\
      Ee=="H"{next}\
      {print}' |\
    cat >! ${t}useme.pdb
else
    cp $pdbfile ${t}useme.pdb
endif

# now we need to mask out the protein
echo "counting conformers in $pdbfile" | tee -a $logfile
cat ${t}useme.pdb |\
awk '! /^ATOM|^HETAT/{next}\
 {conf=substr($0,17,1);occ=substr($0,55,6);\
  if(conf==" ")next;\
  ++count[conf];sum[conf]+=occ}\
 ! order[conf]{++c;order[conf]=c}\
 END{for(conf in count){\
   print order[conf],conf,sum[conf]/count[conf],count[conf]}}' |\
sort -g |\
tee ${t}conf_occ.txt >> $logfile

set masks = ""
foreach cnfnum ( `awk '{print $1}' ${t}conf_occ.txt` )
set cnfchar = `egrep "^${cnfnum} " ${t}conf_occ.txt | awk '{print $2}'`

echo "masking around conf $cnfchar "

egrep "^${cnfnum} " ${t}conf_occ.txt |\
cat - ${t}useme.pdb |\
awk 'NR==1{selected=$2;avgocc=$3+1e-6;next}\
  ! /^ATOM|^HETAT/{print;next}\
    {conf=substr($0,17,1);occ=substr($0,55,6);\
     pre=substr($0,1,54);post=substr($0,61)}\
 conf==" "{print;next}\
 conf!=selected{next}\
    {occ=occ/avgocc}\
  occ>1{occ=1} occ<0{occ=0}\
  {printf("%s%6.2f%s\n",pre,occ,post)}' |\
cat >! maskme_${cnfnum}.pdb
Tc_maskify.com fofc.map maskme_${cnfnum}.pdb outfile=mask_${cnfnum}.map >! ${t}maskify_${cnfnum}.log
set masks = ( $masks mask_${cnfnum}.map )
end

echo "averaging $#masks masks..."
addup_maps_runme.com $masks >&! ${t}addup.log
if( $status || ! -e sum.map) then
  cat ${t}addup.log
  set BAD = "cannot addup_maps"
  goto exit
endif

set scale = `echo $#masks | awk '{print 1/$1}'`

mapmask mapin sum.map mapout not_protein_rough.map << EOF >> $logfile
scale factor $scale
axis $mapAXIS
EOF
#map_func.com -func negate not_protein.map -outfile squish_solvent_rough.map
echo | mapdump mapin not_protein_rough.map | egrep -i "density"

echo "stretching range of protein mask: $protein_lowprob : $protein_highprob -> 0 : 1"
rm -f ${t}temp.map ${t}stretched.map ${t}loclip.map ${t}clipped.map > /dev/null
set mult = `echo $protein_highprob $protein_lowprob | awk '{print 1./($1-$2)}'`
set offs = `echo $protein_highprob $protein_lowprob | awk '{print -$2/($1-$2)}'`
map_func.com -func mult -param $mult not_protein_rough.map -output ${t}temp.map >> $logfile
map_func.com -func add -param $offs ${t}temp.map -output ${t}stretched.map >> $logfile
map_func.com -func max -param 0 ${t}stretched.map -output ${t}loclip.map >> $logfile
map_func.com -func min -param 1 ${t}loclip.map -output ${t}clipped.map >> $logfile

cp ${t}clipped.map not_protein_mask.map
# not protein is a map that is 0  around the protein and 1 far from it
echo | mapdump mapin not_protein_mask.map | egrep -i "density"


echo maps mult |\
mapmask mapin1 not_protein_mask.map \
        mapin2 fofc.map \
        mapout fofc_notprotein.map >> $logfile

skip_protein:





# use 100 bins
set binsize = `echo $reso | awk '{print (0.5/$1)**2/100}'`

#foreach map ( fofc fofc_notprotein Fpart gt_changes FC_ALL_LS mfodfc 2mfodfc )
foreach map ( fofc fofc_notprotein )

if( ! -e "${map}.map" ) continue 

echo "Wilsonifying $map"

gemmi map2sf --dmin=$reso ${map}.map ${t}map.mtz F PHI

set pow = 2
mtzdump hklin ${t}map.mtz << EOF >! ${t}mtzdump.txt
FORMAT '(3i6,3f15.7)'
NREF -1
lreso
go
EOF
# make a Wilson plot
cat ${t}mtzdump.txt |\
awk '/ LIST OF REFLECTIONS/,/ MTZDUMP/' |\
awk 'substr($0,1,18)==sprintf("%6d%6d%6d",$1,$2,$3) && $5+0 !~ -999{print}' |\
tee ${t}mtz.hkl |\
awk -v pow=$pow '{print $4/4,$5**pow}' |\
awk -v bs=$binsize '{bin=sprintf("%.0f",$1/bs);sum1[bin]+=$2;++count[bin]}\
     END{for(bin in sum1) if(bin+0>0) print bin*bs,log(sum1[bin]/count[bin])}' |\
sort -g >! wilson_${map}.txt 
# format: stol^2 log(sumF**2)

#set loline = `awk -v s=$loend '$1<s' wilson_${map}.txt | linfit.awk`
#set hiline = `awk -v s=$hiend '$1>s' wilson_${map}.txt | linfit.awk`

cat << EOF >! fitparams.txt
scale = 0
B = 0
hiscale = 0
hiB = 0
loscale = 0
loB = 0
a5 = 0
a4 = 0
a3 = 0
a2 = 0
a1 = 0
a0 = 0
EOF

gnuplot << EOF >&! ${t}gnuplot.log
reso = $reso
hifrac = $wilsonify_hifrac
lofrac = $wilsonify_lofrac
ssqmax = (0.5/reso)**2
hiend = ssqmax*(1-hifrac)
loend = ssqmax*lofrac
line(x) = scale-2*B*x
hiline(x) = hiscale-2*hiB*x
loline(x) = loscale-2*loB*x
poly(x) = a5*x**5+a4*x**4+a3*x**3+a2*x**2+a1*x+a0
lowt(x) = (x<0?0:(x>loend?1:x/loend))
hiwt(x) = (x>ssqmax?0:(x<hiend?1:1-(x-hiend)/(ssqmax-hiend)))
blend(x) = (1-lowt(x))*loline(x)+(1-hiwt(x))*hiline(x)+lowt(x)*hiwt(x)*poly(x)
# guard: fit the low-res line with a SQUARED slope so loB>=0.  A rising low-res
# Wilson line (loB<0) extrapolates to runaway high-res over-sharpening; this
# happens when an anomalously low origin bin tilts the narrow lofrac fit.
# No-op when the true loB is positive (loBr just converges to sqrt(loB)).
loBr = 5
lolinefit(x) = loscale-2*(loBr*loBr)*x
fit line(x) 'wilson_${map}.txt' via scale,B
fit poly(x) 'wilson_${map}.txt' via a5,a4,a3,a2,a1,a0
fit [:loend] lolinefit(x) 'wilson_${map}.txt' via loscale,loBr
loB = loBr*loBr
fit [hiend:] hiline(x) 'wilson_${map}.txt' via hiscale,hiB
fit blend(x) 'wilson_${map}.txt' via a5,a4,a3,a2,a1,a0
fit blend(x) 'wilson_${map}.txt' via hiscale,loscale,a5,a4,a3,a2,a1,a0
#fit blend(x) 'wilson_${map}.txt' via hiscale,hiB,loscale,loB,a5,a4,a3,a2,a1,a0
set table "${t}table"
set output "${t}table"
set samples 10000
plot [-0.01:ssqmax+0.01] loline(x)-blend(x)
set output
unset table
update "fitparams.txt"
print scale,B,hiscale,hiB,loscale,loB,a5,a4,a3,a2,a1,a0
EOF
set info = `tail -n 1 ${t}gnuplot.log | awk 'BEGIN{RS=" "} {print $0+0}'`
if( $#info != 12) then
    set BAD = "gnuplot fit to Wilson plot failed"
    goto exit
endif
echo "fitparams ${map}: $info"
set newWilsonB = `echo $info | awk '{B=10} $6>B{B=$6} $2>B{B=$2} {print B}'`
echo "highest Wilson B: $newWilsonB"

# ---- build the per-resolution filter (${t}scales.txt lines: "s s s <1/d^2> <H>") ----
set ssqmax = `echo $reso | awk '{print (0.5/$1)**2}'`
if( "$resfilter" == "bandpass" ) then
    echo "resolution filter: band-pass  center(stol^2)=$bandpass_center width=$bandpass_width"
    awk -v s0=$bandpass_center -v ds=$bandpass_width -v smax=$ssqmax \
      'BEGIN{n=4000; for(i=0;i<=n;i++){ss=smax*i/n; H=exp(-((ss-s0)*(ss-s0))/(2.0*ds*ds)); print "s s s",ss*4,H}}' >! ${t}scales.txt
else if( "$resfilter" == "blur" ) then
    echo "resolution filter: blur  B=$blur_B"
    awk -v B=$blur_B -v smax=$ssqmax \
      'BEGIN{n=4000; for(i=0;i<=n;i++){ss=smax*i/n; H=exp(-B*ss); print "s s s",ss*4,H}}' >! ${t}scales.txt
else
    # wilson (default): straighten to the low-res slope, from the gnuplot fit above
    awk '$3=="i" && $2>-38{print "s s s",$1*4,exp($2*0.5)}' ${t}table >! ${t}scales.txt
endif
 sort -k4g ${t}scales.txt ${t}mtz.hkl |\
awk '/^s s s /{scale=$5;next}\
  {print $1,$2,$3,$5*scale,$6,scale}' |\
cat >! ${t}importme.hkl
# format: H K L F phi scale

#echo "re-importing"
f2mtz hklin ${t}importme.hkl hklout ${t}imported.mtz << EOF >> $logfile
SYMM $SGnum
CELL $CELL
LABOUT H K L F PHI
CTYPOUT H H H F P
END
EOF

cad hklin1 ${t}imported.mtz hklout ${t}${map}_Wilsonified.mtz << EOF >> $logfile
labin file 1 all
EOF

fft hklin ${t}${map}_Wilsonified.mtz mapout ${t}ffted.map << EOF >> $logfile
labin F1=F PHI=PHI
$mapGRID
EOF
mapmask mapin ${t}ffted.map mapout ${map}_Wilsonified.map << EOF >> $logfile
xyzlim $mapLIMITS
axis $mapAXIS
EOF

echo "correlate section" |\
overlapmap mapin1 ${map}.map mapin2 ${map}_Wilsonified.map |\
tee ${t}CC.log |\
awk '/Total corr/{print " CC=",$NF}' 


mtzdump hklin ${t}${map}_Wilsonified.mtz << EOF >! ${t}mtzdump.txt
FORMAT '(3i6,3f15.7)'
NREF -1
lreso
go
EOF

# make a new Wilson plot
cat ${t}mtzdump.txt |\
awk '/ LIST OF REFLECTIONS/,/ MTZDUMP/' |\
awk -v pow=$pow 'substr($0,1,18)==sprintf("%6d%6d%6d",$1,$2,$3) && $5+0 !~ -999{\
     print $4/4,$5**pow}' |\
awk -v bs=$binsize '{bin=sprintf("%.0f",$1/bs);sum1[bin]+=$2;++count[bin]}\
     END{for(bin in sum1) if(bin+0>0) print bin*bs,log(sum1[bin]/count[bin])}' |\
sort -g >! wilsonified_${map}.txt 

end


set fofc_style = ""
if( $exclude_protein ) then
    set fofc_style = "${fofc_style}_notprotein"
endif
if( $wilsonify ) then
    set fofc_style = "${fofc_style}_Wilsonified"
endif

if( "$significance" == "none" || "$significance" == "0" ) then
    echo "skipping significance filter"
    set sig_style = ""
    map_func.com -func set -param 1 fofc${fofc_style}.map -outfile significance.map >> $logfile
    goto skip_significance
else
    set sig_style = "significant"
endif


echo "using fofc${fofc_style}.map for significance determination"
echo | mapdump mapin fofc${fofc_style}.map | egrep dens | tee ${t}mapstats.txt
set dens = `awk '/density/{gsub("\\.\\."," ");print $NF}' ${t}mapstats.txt`
set mapstats = `echo $dens | awk '{max=$2} -$1>max{max=-$1} {print max,$4}'`


# automatically determine grow radius from Wilson B and map signal-to-noise ratio
# kernel_d(x) = 4./3*pi/d**3*sinc3(2*pi*x/d)
# sinc3(x) = (x==0?1:3*(sin(x)/x-cos(x))/(x*x))
# kernel_B(x) = (4*pi/B)**1.5*exp(-4*pi**2/B*x*x)

# peak/floor = exp(-4*pi**2*r*r/B)
# log(snr) = 4*pi**2/B*r*r
# log(snr)*B/(4*pi**2) = r*r
# r = sqrt(log(snr)*B/(4.0*pi**2))    
set snr = `echo $mapstats | awk '{print $1/$2}'`
echo "signal/noise = $snr"
echo $newWilsonB $snr 5 1 4.0 |\
awk '{B=$1;snr=$2;snrlim=$3;floor=$4;maxradius=$5;\
      pi=atan2(1,1)*4;}\
      snr>snrlim{snr=snrlim}\
      snr<floor{snr=floor}\
  {radius=sqrt(log(snr)*B/(4.*pi**2))}\
  radius>maxradius{radius=maxradius}\
  {print radius}' >! ${t}autorad.txt
set aradius = `awk '{print $1}' ${t}autorad.txt`

if ( "$radius" == "auto" || $?auto_radius ) then
    set auto_radius
    set radius = $aradius
    if( "$radius" == "" ) set radius = 0
    echo "auto-selected grow radius: $radius"
else
    echo "predicted optimal radius: $aradius"
endif





# ' "probability that something is there" '
# ' prob(rho) = sign(rho)*pow(erf(abs(rho/sigma(rho))/sqrt(2)),V/2/d^3) '




if( "$ShannonB" == "auto" || $?auto_shannon ) then
    set auto_shannon
    echo "auto-setting: ShannonB = $newWilsonB"
    set ShannonB = "$newWilsonB"
endif

# number of Shannon voxels expected
set exponent = `echo $CELLvolume $symops $ShannonB | awk '{V=$1;n=$2;B=$3+$4;pi=atan2(1,1)*4} B<1e-6{B=1e-6} {print V/n/2/(sqrt(B*log(2))/pi)**3}'`
echo "estimating $exponent Shannon voxels in ASU"
if( "$exponent" == "" || "$exponent" == "0" ) then
    set BAD = "error computing number of Shannon voxels in map"
    goto exit
endif

#set goal_sigma = "0.7071"
set maxprob = `echo "print erf(abs(${snr}*sqrt(2)))**$exponent" | gnuplot |& cat`
set minsnr = `echo "print inverf(($highprob-1e-6)**(1./$exponent))/sqrt(2)" | gnuplot |& cat`
set goal_sigma = `echo $snr $minsnr | awk '{su=$2/$1;gs=sqrt(0.5)} su<1{su=1} {print su*gs}'`
echo "goal sigma: $goal_sigma  maxprob: $maxprob   minsnr: $minsnr"

echo "probability unmodeled = erf(abs(rho/sigma* sqrt(2)))^$exponent : $maxprob"
mapmask mapin fofc${fofc_style}.map mapout ${t}fofc_sigscaled.map << EOF >> $logfile
scale sigma $goal_sigma 0
EOF
map_func.com -func abs  ${t}fofc_sigscaled.map -output ${t}abs.map >> $logfile
map_func.com -func erfpow -param $exponent  ${t}abs.map -output ${t}erfpow.map >> $logfile
if( $status ) then
    set BAD = "error computing erfpow of map"
    goto exit
endif
cp ${t}erfpow.map rough_significance.map
echo "sharp-edged significance mask in: rough_significance.map"
echo | mapdump mapin rough_significance.map | egrep density

echo "padding map..."
echo | mapdump mapin ${t}erfpow.map >! ${t}mapdump.txt
set xyzlim = `awk '/Fast, medium, slow/{o[$(NF-2)]=1;o[$(NF-1)]=2;o[$NF]=3; print b[o["X"]],e[o["X"]], b[o["Y"]],e[o["Y"]], b[o["Z"]],e[o["Z"]]; exit} /Start and stop points/{b[1]=$(NF-5); e[1]=$(NF-4); b[2]=$(NF-3); e[2]=$(NF-2); b[3]=$(NF-1); e[3]=$NF}' ${t}mapdump.txt`

set radius_pixel = `echo $radius $voxel_size | awk '{print int($1/$NF)}'`
set pixrad3d = `echo $radius $voxel_size | awk '{print "-xrad",$1/$2,"-yrad",$1/$3,"-zrad",$1/$4}'`
set radius_residual = `echo $radius $radius_pixel $voxel_size | awk '{print $1-$2*$NF}'`
if( $debug ) echo "DEBUG: radius_pixel: $radius_pixel  resid: $radius_residual"

set bevel_pixel = `echo $bevel $voxel_size | awk '{print int($1/$NF)}'`
set pixbevel3d = `echo $bevel $voxel_size | awk '{print "-xrad",$1/$2,"-yrad",$1/$3,"-zrad",$1/$4}'`
set bevel_residual = `echo $bevel $bevel_pixel $voxel_size | awk '{print $1-$2*$NF}'`
if( $debug ) echo "DEBUG: bevel_pixel: $bevel_pixel  resid: $bevel_residual"
set padding = `echo $radius_pixel $bevel_pixel 5 | awk '{print $1+$2+$3}'`
if( $debug ) echo "DEBUG: padding by $padding voxels in all directions"


set paddedlim = `echo $xyzlim $padding | awk '{p=$NF;print $1-p,$2+p,$3-p,$4+p,$5-p,$6+p}'`
if( $debug ) echo "DEBUG: xyzlim: $xyzlim"
if( $debug ) echo "DEBUG: xyzlim: $paddedlim"

mapmask mapin ${t}erfpow.map mapout ${t}padded.map << EOF >> $logfile
xyzlim $paddedlim
EOF
if( $status ) then
    set BAD = "error extending map for preffting"
    goto exit
endif
echo "growing selection by $radius A ($radius_pixel voxels)"

set grow_rounds = `echo $radius_pixel $max_pix_grow | awk '{print int($1/$2)+1}'`
#set grow_rounds = 1
set incr_radius = `echo $max_pix_grow $voxel_size | awk '{print $1*$2}'`
set final_radius = `echo $radius $grow_rounds $incr_radius | awk '{print $1-($2-1)*$3}'`
if( $debug ) echo "DEBUG: grow_rounds = $grow_rounds incr_radius = $incr_radius final_radius = $final_radius "
foreach mr ( `seq 1 $grow_rounds` )
  if( $mr == $grow_rounds ) set incr_radius = $final_radius
  set pixrad3d = `echo $incr_radius $voxel_size | awk '{print "-xrad",$1/$2,"-yrad",$1/$3,"-zrad",$1/$4}'`
  echo "growing selection by $incr_radius A ($pixrad3d[2] voxels) this round ( $mr of $grow_rounds )"
  map_func.com -func maxradius $pixrad3d ${t}padded.map output=${t}padmaxrad.map  >> $logfile
  cp ${t}padmaxrad.map ${t}padded.map
end
map_func.com -func maxradius $pixbevel3d ${t}padded.map output=${t}padbeveled.map  >> $logfile
# always add an extra 2 voxels for "zero" mask
map_func.com -func maxradius 2 ${t}padbeveled.map output=${t}padmask.map  >> $logfile

# save thresholded versions for masking after fft
mapmask mapin ${t}padmask.map mapout ${t}minmask.map << EOF >> $logfile
xyzlim $mapLIMITS
EOF
mapmask mapin ${t}padmaxrad.map mapout ${t}maxrad.map << EOF >> $logfile
xyzlim $mapLIMITS
EOF
map_func.com -func thresh -param $highprob ${t}erfpow.map -output ${t}stayup.map >> $logfile
map_func.com -func thresh -param $lowprob ${t}minmask.map -output ${t}notlow.map >> $logfile



echo "pre-fft conditioning of sharp edges"
map_func.com -func prefft ${t}padmaxrad.map -output ${t}padprefft.map >> $logfile
mapmask mapin ${t}padprefft.map mapout ${t}fftme.map << EOF >> $logfile
xyzlim $mapLIMITS
EOF
cp ${t}fftme.map preffted_significance.map
echo | mapdump mapin preffted_significance.map | egrep density

foreach itr ( `seq 1 $recycles` )
    set thisB = $filterB
    set descr = preffted
    if( $itr > 1 ) set descr = reinforced
    if( $itr == $recycles ) then
      set thisB = 1e-6
    endif
    echo -n "recombining $descr mask with ffted version of itself: $itr of $recycles"
    map_scaleB.com ${t}fftme.map B=$thisB output=${t}filtered.map clip=0 >> $logfile
    if( $status ) then
        set BAD = "error re-ffting map"
        goto exit
    endif

    # make sure selected areas stay selected
    map_func.com -func max ${t}stayup.map ${t}filtered.map -output ${t}maxmix.map >> $logfile
    # make sure de-selected areas stay de-selected
    map_func.com -func mult ${t}notlow.map ${t}maxmix.map -output ${t}reinforced.map >> $logfile
#    cp ${t}maxmix.map ${t}reinforced.map

    # do some range stretching to get rid of ripples
    set mult = `echo 1 0 $highprob $lowprob | awk '{print ($1-$2)/($3-$4)}'`
    set offs = `echo $mult 0 $lowprob | awk '{print ($2-$3)*$1}'`
    if( $debug > 3 ) echo "\nDEBUG: mult= $mult  offs= $offs"
    map_func.com -func mult -param $mult ${t}reinforced.map -output ${t}mult.map >> $logfile
    map_func.com -func add -param $offs ${t}mult.map -output ${t}stretched.map >> $logfile
    map_func.com -func max -param 0 ${t}stretched.map -output ${t}loclip.map >> $logfile
    map_func.com -func min -param 1 ${t}loclip.map -output ${t}clipped.map >> $logfile


    echo "correlate section" |\
    overlapmap mapin1 ${t}filtered.map mapin2 ${t}fftme.map |\
    tee ${t}CC.log |\
    awk '/Total corr/{print " CC=",$NF}' 

    cp ${t}clipped.map ${t}fftme.map
end
cp ${t}clipped.map filtered_significance.map
echo "filtered_significance.map will now fft with minimal aliasing noise"
echo | mapdump mapin filtered_significance.map | egrep density

# do something here to prevent radius and bevel from getting too small
# or from radius from getting too large relative to bevel

# compute scale and B needed to generate specified radius and bevel
echo "beveling by $bevel A" | tee -a $logfile

# pick a B factor that has fwhm = bevel
echo $bevel |\
awk '{bevel=$1;pi=4*atan2(1,1);\
  B = -4*pi**2*(bevel/2)**2/log(0.5);\
  s=(4*pi/B)**1.5;\
  print s,B}' >! ${t}sB.txt
set bevel_scale = `awk 'NF==2{print $1}' ${t}sB.txt`
set bevel_B     = `awk 'NF==2{print $2}' ${t}sB.txt`
if( $debug ) echo "DEBUG: bevel_scale= $bevel_scale bevel_B= $bevel_B"
if( "$bevel_B" == "") then
    set BAD = "unable to compute bevel B"
    goto exit
endif

cp filtered_significance.map ${t}fftme.map

foreach itr ( `seq 1 $recycles` )
    set thisB = $bevel_B
    set descr = beveled
    if( $itr > 1 ) set descr = reinforced
    if( $itr == $recycles ) then
      set thisB = 1e-6
    endif
    echo -n "recombining $descr mask with ffted version of itself: $itr of $recycles"
    map_scaleB.com ${t}fftme.map B=$thisB output=${t}filtered.map >> $logfile
    if( $status ) then
        set BAD = "error re-ffting map"
        goto exit
    endif

    # make sure selected areas stay selected
    map_func.com -func max ${t}stayup.map ${t}filtered.map -output ${t}maxmix.map >> $logfile
    # make sure de-selected areas stay de-selected
    map_func.com -func mult ${t}notlow.map ${t}maxmix.map -output ${t}reinforced.map >> $logfile

    # do some range stretching to get rid of ripples
    set mult = `echo 1 0 $highprob $lowprob | awk '{print ($1-$2)/($3-$4)}'`
    set offs = `echo $mult 0 $lowprob | awk '{print ($2-$3)*$1}'`
    if( $debug ) echo "\nDEBUG: mult= $mult  offs= $offs"
    map_func.com -func mult -param $mult ${t}reinforced.map -output ${t}mult.map >> $logfile
    map_func.com -func add -param $offs ${t}mult.map -output ${t}stretched.map >> $logfile
    map_func.com -func max -param 0 ${t}stretched.map -output ${t}loclip.map >> $logfile
    map_func.com -func min -param 1 ${t}loclip.map -output ${t}clipped.map >> $logfile

    echo "correlate section" |\
    overlapmap mapin1 ${t}filtered.map mapin2 ${t}fftme.map |\
    tee ${t}CC.log |\
    awk '/Total corr/{print " CC=",$NF}' 

    cp ${t}clipped.map ${t}fftme.map
end

cp ${t}clipped.map significance.map

echo "significance map reflects probability that the fo-fc feature is not noise"
echo | mapdump mapin significance.map | egrep density




# now analyze before-and-after histograms
set n = 0
foreach map ( rough_significance filtered_significance significance )
    @ n = ( $n + 1 )
    if( $debug ) echo "DEBUG: hist$n = $map"
    dd bs=$mapheadersize skip=1 if=${map}.map of=map.bin >& /dev/null
    float_add map.bin -histogram -binsize 0.01 -outfile /dev/null >! ${t}hist.log
    awk '/^histogram/,$2=="pixels"' ${t}hist.log | awk '! /[a-z]/' >! ${t}_hist${n}.txt
    cp ${t}_hist${n}.txt ${map}_hist.txt
end
cp ${t}_hist1.txt ${t}before_hist.txt
cp ${t}_hist${n}.txt ${t}after_hist.txt

if( ! -s ${t}_hist${n}.txt) then
  # try using iotbx
  iotbx.python << EOF >! ${t}hists.txt
from iotbx import ccp4_map
import numpy as np

m = ccp4_map.map_reader("filtered_significance.map")
vals = m.data.as_1d()
hist, bins = np.histogram(vals, bins=100)
for i in range(100):
    print(f"{bins[i]:8.5g} {hist[i]:8d}")
m = ccp4_map.map_reader("significance.map")
vals = m.data.as_1d()
hist, bins = np.histogram(vals, bins=100)
for i in range(100):
    print(f"{bins[i]:8.5g} {hist[i]:8d}")
EOF
  head -n 100 ${t}hists.txt >! ${t}before_hist.txt
  tail -n 100 ${t}hists.txt >! ${t}after_hist.txt
endif
set locut0 = `awk '{rho[NR]=$1;N[NR]=$2} $2>max{maxv=$1;max=$2} $2<max{print maxv;exit}' ${t}before_hist.txt`
set hicut0 = `tac ${t}before_hist.txt | awk '{rho[NR]=$1;N[NR]=$2} $2>max{maxv=$1;max=$2} $2<max{print maxv;exit}'`
awk '$1<0.7 && $2>1000' ${t}after_hist.txt |\
awk '{rho[NR]=$1;count[NR]=$2;prev=NR-1;prevprev=NR-2}\
  count[prev]>count[prevprev] && count[prev]>count[NR]{\
   max=count[prev];\
   print rho[prev],max,max/((count[prevprev]+count[NR])/2)}' >! ${t}peaks.txt
set locut = `awk '$3>1.1{print $1}' ${t}peaks.txt | head -n 1`
tac ${t}after_hist.txt |\
awk '$1>0.2 && $2>1000' |\
awk '{rho[NR]=$1;count[NR]=$2;prev=NR-1;prevprev=NR-2}\
  count[prev]>count[prevprev] && count[prev]>count[NR]{\
   max=count[prev];\
   print rho[prev],max,max/((count[prevprev]+count[NR])/2)}' >! ${t}peaks.txt
set hicut = `head -1 ${t}peaks.txt | awk '{print $1}'`
if( $debug ) echo "DEBUG: cuts: $locut0 $hicut0  $locut $hicut"
if( "$locut" == "" ) set locut = 0
if( "$hicut" == "" ) set hicut = 1
echo "recommend stretching beveled map from current range ( $locut : $hicut ) to match pre-blur range ( $locut0 : $hicut0 )"


skip_significance:



if( ! $exclude_protein ) then
    echo "not excluding protein region"
    cp significance.map significant_solvent_mask.map 
else
    echo "merging mask of significant difference density with the is-not-protein mask"
    echo maps mult |\
    mapmask mapin1 not_protein_mask.map \
            mapin2 significance.map \
            mapout significant_solvent_mask.map >> $logfile
endif
# mask of solvent differences that are significant




if( "$nearpeak" == "none" || "$nearpeak" == "0" ) then
    if( $exclude_protein ) set sig_style = "$sig_style non-protein"
    echo "using all ${sig_style} features"
    goto skip_nearpeak
endif

# need to make signif map so we can pick a peak
echo "fofc_signif_sol is significant solvent features only"
map_func.com -func multiply fofc${fofc_style}.map significant_solvent_mask.map -outfile fofc_signif_sol_full.map >> $logfile
echo | mapdump mapin fofc_signif_sol_full.map | egrep -i "density"
 
# compute scale and B needed to generate specified radius and bevel
echo "selecting area within $nearpeak_radius A of maximum, beveled by $nearpeak_bevel A"
echo $nearpeak_radius $nearpeak_bevel |\
awk '{radius=$1;bevel=$2;pi=4*atan2(1,1);\
  rat = ((radius+bevel/2)/radius)**2;\
  s = (0.977322 - 0.159175*rat)/(rat - 1);\
  B = -4*pi**2*radius**2/log(-log(0.5)*10**-s);\
  print -(10**s),B}' >! ${t}sB.txt
set bevel_scale = `awk '{print $1}' ${t}sB.txt`
set bevel_B     = `awk '{print $2}' ${t}sB.txt`
if("$bevel_scale" == "" || "$bevel_B" == "") then
    set BAD = "radius and bevel combination too sharp"
    goto exit
endif
if( $debug ) echo "DEBUG: bevel_scale= $bevel_scale bevel_B= $bevel_B"

# calculate B factor that will already be applied by resolution cutoff
set resoB = `echo $reso | awk '{print 9*$1*$1}'`


set pickatom_B = `echo $bevel_B | awk '$1>100{$1=100} {print $1}'`
set residual_B = `echo $bevel_B $pickatom_B | awk '{print $1-$2}'`

# highest point in map
echo | mapdump mapin fofc_signif_sol_full.map >! ${t}mapdump.txt
set stats = `awk '/density/{gsub("\\.\\."," ");print $NF}' ${t}mapdump.txt`
if( $debug ) echo "DEBUG: fofc_signif_sol density stats: $stats"

set top_bottom = `echo $stats | awk '$2>=-$1{print "top";exit} $2<-$1{print "bottom"}' | tail -n 1`
if( $debug ) echo "DEBUG: top_botttom = $top_bottom "
pick.com -$top_bottom fofc_signif_sol_full.map  >> $logfile
if( $status ) then
    set BAD = "peak pick failed"
    goto exit
endif
echo $pickatom_B |\
cat - pick.pdb |\
awk 'NR==1{B=$1;next} ! /^ATOM/{print;next} \
  {printf("ATOM%7d ANO  ANO A   1    %s  1.00%6.2f\n",1, substr($0,31,24),B)}' |\
cat >! ${t}pickatom.pdb
awk '/^ATOM/{print "peak density at:",substr($0,31,24)}' ${t}pickatom.pdb
set test = `egrep "^ATOM|^HETAT" ${t}pickatom.pdb | wc -l`
if( "$test" == "0" || "$test" == "" ) then
    set BAD = "no peaks found"
    goto exit
endif

sfall xyzin ${t}pickatom.pdb mapout ${t}sfalled.map << EOF | tee ${t}sfall.log  >> $logfile
mode atmmap
SYMM $SGnum
$mapGRID
#vdwr 3
SFSG 1
EOF
#gemmi sfcalc -v  --dmin=$reso --rate 2.8 \
#     --write-map=${t}gemmi.map --to-mtz=${t}gemmi.mtz ${t}pickatom.pdb

mapmask mapin ${t}sfalled.map mapout ${t}raw.map << EOF >> $logfile
xyzlim $mapLIMITS
axis $mapAXIS
EOF
map_scaleB_runme.com ${t}raw.map output=${t}bigB.map \
    B=$residual_B reso=$reso clip=0  >> $logfile
echo "scale factor $bevel_scale" |\
mapmask mapin ${t}bigB.map mapout ${t}neg.map >> $logfile
map_func.com -func min -param 0 ${t}neg.map -output ${t}expme.map >> $logfile
map_func.com -func safexp ${t}expme.map -output ${t}exp.map >> $logfile
map_func.com -func negate -param 1 ${t}exp.map -output nearpeak_mask_rough.map >> $logfile
cp nearpeak_mask_rough.map nearpeak_mask.map

foreach itr ( `seq 1 $nearpeak_recycles` )
    echo -n "recombining near-peak mask with ffted version of itself: $itr of $nearpeak_recycles"
    map_scaleB.com nearpeak_mask.map B=$filterB output=${t}filtered.map clip=0 >> $logfile
    if( $status ) then
        set BAD = "error re-ffting map"
        goto exit
    endif

    map_func.com -func max nearpeak_mask.map ${t}filtered.map -output ${t}maxmix.map >> $logfile
    map_func.com -func min -param 1 ${t}maxmix.map -output ${t}clipped.map >> $logfile

    echo "correlate section" |\
    overlapmap mapin1 ${t}clipped.map mapin2 nearpeak_mask.map |\
    tee ${t}CC.log |\
    awk '/Total corr/{print " CC=",$NF}' 

    mv ${t}clipped.map nearpeak_mask.map
end



echo "merging nearkpeak_mask with significant difference mask"
mv significant_solvent_mask.map full_significant_solvent_mask.map
echo maps mult |\
mapmask mapin1 nearpeak_mask.map mapin2 full_significant_solvent_mask.map mapout significant_solvent_mask.map >> $logfile

skip_nearpeak:



# implement segmentation-based flood-filling around peak
if( "$floodfill" == "none" || "$floodfill" == "0" ) then
    echo "skipping floodfill"
    goto skip_floodfill
endif

echo "flood-filling around most-extreme peak in fofc${fofc_style}.map"
if(! -e fofc${fofc_style}.map) then
  set BAD = "missed a step"
  goto exit
endif

# highest/lowest point in map
echo | mapdump mapin fofc${fofc_style}.map >! ${t}mapdump.txt
set stats = `awk '/density/{gsub("\\.\\."," ");print $NF}' ${t}mapdump.txt`
if( $debug ) echo "DEBUG: fofc_signif_sol density stats: $stats"

set top_bottom = `echo $stats | awk '$2>=-$1{print "top";exit} $2<-$1{print "bottom"}' | tail -n 1`
if( $debug ) echo "DEBUG: top_botttom = $top_bottom "
pick.com -$top_bottom fofc${fofc_style}.map  >> $logfile
if( $status ) then
    set BAD = "peak pick failed"
    goto exit
endif
cat pick.pdb |\
awk '! /^ATOM/{print;next} \
  {printf("ATOM%7d ANO  ANO A   1    %s  1.00%6.2f\n",1, substr($0,31,24),B)}' |\
cat >! ${t}pickatom.pdb
awk '/^ATOM/{print "peak density at:",substr($0,31,24)}' ${t}pickatom.pdb
set test = `egrep "^ATOM|^HETAT" ${t}pickatom.pdb | wc -l`
if( "$test" == "0" || "$test" == "" ) then
    set BAD = "no peaks found"
    goto exit
endif

cp fofc${fofc_style}.map ${t}smoothme.map
if( "$top_bottom" == "bottom" ) then
   echo scale factor -1 |\
   mapmask mapin1 fofc${fofc_style}.map mapout ${t}smoothme.map >> $logfile
endif

set edgevoxels = 1
while ( $edgevoxels )

echo -n "settings: "
cat << EOF | awk 'BEGIN{ORS=" "} {print}'
floodfill_sigma=$floodfill_sigma
floodfill_Bsmooth=$floodfill_Bsmooth
floodfill_maskBsmooth=$floodfill_maskBsmooth
floodfill_rethresh=$floodfill_rethresh
floodfill_maxradius=$floodfill_maxradius
floodfill_recycles=$floodfill_recycles
EOF
echo ""

echo "smoothing fofc${fofc_style} with floodfill_Bsmooth=$floodfill_Bsmooth"
map_scaleB.com B=$floodfill_Bsmooth \
  ${t}smoothme.map outfile=${t}threshme.map >> $logfile
echo | mapdump mapin ${t}threshme.map | egrep density >! ${t}smoothstat.txt

echo | mapdump mapin ${t}threshme.map >! ${t}mapdump.txt
set smooth_stats = `awk '/density/{print $NF}' ${t}mapdump.txt`
set smooth_rmsd = `echo $smooth_stats | awk '{print $4}'`

set pad_radius = `echo $floodfill_maxradius 1 | awk '{print $1+$2}'`
echo "padding edges of flood-fill box by 1 A : $pad_radius"
mapmask mapin1 ${t}threshme.map xyzin pick.pdb mapout ${t}chunk.map << EOF >> $logfile
border $pad_radius
EOF

# initial threshold - for non-mask map
set thresh = `echo $smooth_rmsd $floodfill_sigma | awk '{print $1*$2}'`
echo "thresholding smooth fofc map at $thresh ($floodfill_sigma sigma)"
map_func.com -func thresh -param $thresh ${t}chunk.map -output ${t}thresh.map >> $logfile

foreach flooditr ( `seq 1 $floodfill_recycles` )

  cp ${t}thresh.map ${t}thresh${flooditr}.map

  # actual segmentation run - requires debugged float_func.c
  echo "segmenting $pad_radius A box  around peak (pass $flooditr of $floodfill_recycles)"
  map_func.com -func segment \
    ${t}thresh.map -output ${t}segments${flooditr}.map | tee ${t}flood_seg${flooditr}.log >> $logfile
  set nsegments = `awk '/Maximum density/{print $NF+0}' ${t}flood_seg${flooditr}.log | tail -n 1`
  if( $debug ) then
    echo "found $nsegments segments"
    # map2pdb.com  ${t}segments${flooditr}.map 
  endif

  # just get a few voxels around the picked point
  set small  = `echo $voxel_size | awk '{print $NF}'`
  echo border $small |\
   mapmask mapin1 ${t}segments${flooditr}.map xyzin pick.pdb mapout ${t}peek.map >> $logfile

  # extract the segment number of the central
  od -f -v -w4 --skip $mapheadersize ${t}peek.map |\
  awk 'NF==2 && $NF+0>0{print $2}' |\
  sort -g |\
   awk '{v[NR]=$1} END{print v[int(NR/2)]}' |\
  cat >! ${t}segnumber.txt
  set segnum = `cat ${t}segnumber.txt`
  echo "peak is in segment number $segnum of $nsegments"
  if( "$segnum" == "" || "$segnum" == "0" ) then
    set BAD = "peak was not in a segment"
    goto exit
  endif

  # isolate segment by subtracting its value and removing everything non-zero
  if( $debug ) echo "isolating segment $segnum"
  map_func.com -func subtract -param $segnum ${t}segments${flooditr}.map -output ${t}zeroseg.map >> $logfile
  map_func.com -func abs ${t}zeroseg.map -output ${t}notseg.map >> $logfile
  map_func.com -func negate  ${t}notseg.map -output ${t}hiseg.map >> $logfile
  map_func.com -func thresh -param 0.5 ${t}hiseg.map -output ${t}thresh_segment.map >> $logfile



  # now do distance transform and watershed segmentation
  # - requires latest float_func.c
  echo "distrans-watershed segmenting ..."
  map_func.com -func disttrans \
    ${t}thresh_segment.map -output ${t}disttrans.map  >> $logfile
  map_func.com -func watershed \
    ${t}disttrans.map -output ${t}watershed.map | tee ${t}watershed_seg${flooditr}.log >> $logfile
  set nsegments = `awk '/Maximum density/{print $NF+0}' ${t}watershed_seg${flooditr}.log | tail -n 1`
  if( $debug ) then
    echo "found $nsegments segments"
    # map2pdb.com  ${t}segments${flooditr}.map 
  endif

  # just get a few voxels around the picked point
  set small  = `echo $voxel_size | awk '{print $NF}'`
  echo border $small |\
   mapmask mapin1 ${t}watershed.map xyzin pick.pdb mapout ${t}peek.map >> $logfile

  # extract the segment number of the central
  od -f -v -w4 --skip $mapheadersize ${t}peek.map |\
  awk 'NF==2 && $NF+0>0{print $2}' |\
  sort -g |\
   awk '{v[NR]=$1} END{print v[int(NR/2)]}' |\
  cat >! ${t}segnumber.txt
  set segnum = `cat ${t}segnumber.txt`
  echo "peak is in segment number $segnum of $nsegments"
  if( "$segnum" == "" || "$segnum" == "0" ) then
    set BAD = "peak was not in a segment"
    goto exit
  endif

  # isolate segment by subtracting its value and removing everything non-zero
  if( $debug ) echo "isolating segment $segnum"
  map_func.com -func subtract -param $segnum ${t}watershed.map -output ${t}zeroseg.map >> $logfile
  map_func.com -func abs ${t}zeroseg.map -output ${t}notseg.map >> $logfile
  map_func.com -func negate  ${t}notseg.map -output ${t}hiseg.map >> $logfile
  map_func.com -func thresh -param 0.5 ${t}hiseg.map -output ${t}segment.map >> $logfile




  # put chunk back into an empty cell
  mapmask mapin ${t}segment.map mapout ${t}extended.map << EOF >> $logfile
  xyzlim $mapLIMITS
  pad 0
EOF


  # run it through FFT to smooth it - helps isolate edges
  # but first, must condition mask for surviving fft
  cp ${t}extended.map sharpedged_floodfill.map

  map_func.com -func prefft sharpedged_floodfill.map -output ${t}fftme.map >> $logfile
  if( $debug ) echo | mapdump mapin ${t}fftme.map | egrep density

  # establish regions to keep re-setting back to 1 and 0
  echo "floodfill fft margin between stayup and not-low masks: $floodfill_fftmargin voxels"
  cp sharpedged_floodfill.map ${t}stayup.map
  map_func.com -func maxradius -param $floodfill_fftmargin ${t}stayup.map -output ${t}notlow.map >> $logfile


  foreach itr ( `seq 1 $recycles` )
      set thisB = `echo $floodfill_filterB $itr | awk '{print $1*(0.9**$2)}'`
      set descr = preffted
      if( $itr > 1 ) set descr = reinforced
      if( $itr == $recycles ) then
        set thisB = 1e-6
      endif
      echo -n "recombining $descr mask with ffted version of itself: $itr of $recycles (B=$thisB)"
      map_scaleB.com ${t}fftme.map B=$thisB output=${t}filtered.map clip=0 >> $logfile
      if( $status ) then
          set BAD = "error re-ffting map"
          goto exit
      endif

      # make sure selected areas stay selected
      map_func.com -func max ${t}stayup.map ${t}filtered.map -output ${t}maxmix.map >> $logfile
      # make sure most de-selected areas stay de-selected
      map_func.com -func mult ${t}notlow.map ${t}maxmix.map -output ${t}reinforced.map >> $logfile

      # do some range stretching to get rid of ripples
      set mult = `echo 1 0 $floodfill_highprob $floodfill_lowprob | awk '{print ($1-$2)/($3-$4)}'`
      set offs = `echo $mult 0 $floodfill_lowprob | awk '{print ($2-$3)*$1}'`
      if( $debug > 1 ) echo "\nDEBUG: mult= $mult  offs= $offs"
      map_func.com -func mult -param $mult ${t}reinforced.map -output ${t}mult.map >> $logfile
      map_func.com -func add -param $offs ${t}mult.map -output ${t}stretched.map >> $logfile
      map_func.com -func max -param 0 ${t}stretched.map -output ${t}loclip.map >> $logfile
      map_func.com -func min -param 1 ${t}loclip.map -output ${t}clipped.map >> $logfile

      echo "correlate section" |\
      overlapmap mapin1 ${t}filtered.map mapin2 ${t}fftme.map |\
      tee ${t}CC.log |\
      awk '/Total corr/{print " CC=",$NF}' 

      cp ${t}clipped.map ${t}fftme.map
  end

  cp ${t}clipped.map filtered_floodfill.map
  echo "filtered_floodfill.map will now fft with minimal aliasing noise"
  if( $debug ) echo | mapdump mapin filtered_floodfill.map | egrep density

  echo border $small |\
   mapmask mapin1 filtered_floodfill.map xyzin pick.pdb mapout ${t}peek.map >> $logfile
  echo | mapdump mapin ${t}peek.map | egrep density >! ${t}before.txt
  set peakpeek = `awk '/Mean dens/{print $NF} /Rms devia/{print $NF}' ${t}before.txt`
  if( $debug ) echo "peak peek: $peakpeek"

  # run it through FFT to smooth it - helps isolate edges
  echo "smoothing mask with floodfill_maskBsmooth=$floodfill_maskBsmooth"
  map_scaleB.com B=$floodfill_Bsmooth filtered_floodfill.map outfile=${t}smoothed.map >> $logfile
  echo | mapdump mapin ${t}smoothed.map | egrep density >! ${t}smoothstat.txt

  echo border $small |\
   mapmask mapin1 ${t}smoothed.map xyzin pick.pdb mapout ${t}peek.map >> $logfile
  echo | mapdump mapin ${t}peek.map | egrep density >! ${t}after.txt
  set peakpeek = `awk '/Mean dens/{print $NF} /Rms devia/{print $NF}' ${t}after.txt`
  if( $debug ) echo "peak peek: $peakpeek"


  # preserve peak area level
  set before = `awk '/Mean dens/{mean=$NF} /Rms devia/{print $NF+mean}' ${t}before.txt`
  set after  = `awk '/Mean dens/{mean=$NF} /Rms devia/{print mean-$NF}' ${t}after.txt`
  set scale = `echo $before $after | awk '$2+0>0{print $1/$2}'`
  if( "$scale" == "" ) then
    set BAD = "lost too many peak voxels in floodfill"
    goto exit
  endif

  echo "scaling by $scale to match peak-area levels"
  echo scale factor $scale |\
  mapmask mapin ${t}smoothed.map mapout ${t}scaled.map >> $logfile
  if( $debug ) echo | mapdump mapin ${t}scaled.map | egrep density

  # do some range stretching to get rid of ripples
  echo "stretching range again"
  set mult = `echo 1 0 $floodfill_highprob $floodfill_lowprob | awk '{print ($1-$2)/($3-$4)}'`
  set offs = `echo $mult 0 $floodfill_lowprob | awk '{print ($2-$3)*$1}'`
  if( $debug > 1 ) echo "\nDEBUG: mult= $mult  offs= $offs"
  map_func.com -func mult -param $mult ${t}scaled.map -output ${t}mult.map >> $logfile
  map_func.com -func add -param $offs ${t}mult.map -output ${t}stretched.map >> $logfile
  map_func.com -func max -param 0 ${t}stretched.map -output ${t}loclip.map >> $logfile
  map_func.com -func min -param 1 ${t}loclip.map -output ${t}clipped2.map >> $logfile
  if( $debug ) echo | mapdump mapin ${t}clipped2.map | egrep density

  # select box around the peak again
  echo "re-defining mask with threshold floodfill_rethresh=$floodfill_rethresh"
  echo border $pad_radius |\
  mapmask mapin1 ${t}clipped2.map xyzin pick.pdb mapout ${t}chunk.map >> $logfile
  map_func.com -func thresh -param $floodfill_rethresh \
     ${t}chunk.map -output ${t}thresh.map >> $logfile

  # check that peak is still selected
  echo border $small |\
   mapmask mapin1 ${t}clipped2.map xyzin pick.pdb mapout ${t}peek.map >> $logfile
  echo | mapdump mapin ${t}peek.map | grep dens >! ${t}peakpeek.txt
  set peakpeek = `awk '/Mean dens/{print $NF} /Rms devia/{print $NF}' ${t}peakpeek.txt`
  if( $debug ) echo "peak peek: $peakpeek"
  set test = `echo $peakpeek | awk '{print ( $1>=0.99 && $2<$1/5 )}'`
  if( ! $test ) then
    set BAD = "peak selection too rough"
    goto exit
  endif

end

# check for edges
mapmask mapin ${t}thresh.map xyzin pick.pdb mapout ${t}dumpme.map << EOF | tee ${t}mapdump.log | grep dens
axis X Y Z
border $floodfill_maxradius
EOF
cat ${t}mapdump.log |\
awk '/Grid sampling on x, y, z/{gx=$8;gy=$9;gz=$10}\
     /Cell dimensions /{xs=$4/gx;ys=$5/gy;zs=$6/gz}\
     /Maximum density /{max=$NF}\
     /Minimum density /{min=$NF}\
     /Start and stop points on cols/{cs=$9+0;rs=$11+0;ss=$13+0}\
     /Number of columns, rows, sections/{nc=$7;nr=$8;ns=$9}\
 END{print xs,ys,zs,cs,rs,ss,nc,nr,ns,gx,gy,gz,max,min}' >! ${t}mapstuff.txt

od -vf -w4 -j $mapheadersize ${tempfile}dumpme.map | awk 'NF==2{print $2}'|\
cat ${tempfile}mapstuff.txt - |\
awk 'NR==1{cs=$4;rs=$5;ss=$6;nc=$7;nr=$8;ns=$9;gx=$10;gy=$11;gz=$12;max=$13;min=$14;\
   c=r=s=0;\
   next}\
  $1>0 && (c==0 || r==0 || s==0 || c==nc-1 || r==nr-1 || s==ns-1){print c,r,s,$1}\
  {++c} c>=nc{c=0;++r} r>=nr{r=0;++s}' |\
cat >! ${t}edges.txt 
set edgevoxels = `cat ${t}edges.txt | wc -l`
echo "$edgevoxels voxels on edge of $floodfill_maxradius A box still selected"

if( ! $edgevoxels ) then
  echo "peak successfully isolated "
  break
endif

if( $floodfill_retries ) then
    @ floodfill_retries = ( $floodfill_retries - 1 )
endif
if( ! $floodfill_retries ) then
    set BAD = "flood-fill hit edge of box"
    goto exit
endif
echo "trying again..."

set test = `echo $floodfill_fftmargin | awk '{print ( $1 < 5 )}'`
if( $test ) then
    set floodfill_fftmargin = `echo $floodfill_fftmargin | awk '{print $1+1}'`
    echo "increasing to floodfill_fftmargin=$floodfill_fftmargin"
    continue
endif

set test = `echo $floodfill_maskBsmooth | awk '{print ( $1 < 20 )}'`
if( $test ) then
    set floodfill_maskBsmooth = `echo $floodfill_maskBsmooth | awk '{print $1*2}'`
    echo "increasing to floodfill_maskBsmooth=$floodfill_maskBsmooth"
    continue
endif

set test = `echo $floodfill_lowprob | awk '{print ( $1 < 0.5 )}'`
if( $test ) then
    set floodfill_lowprob = `echo $floodfill_lowprob | awk '{print $1+0.1}'`
    echo "increasing to floodfill_lowprob=$floodfill_lowprob"
    continue
endif

set test = `echo $floodfill_highprob | awk '{print ( $1 > 0.6 )}'`
if( $test ) then
    set floodfill_highprob = `echo $floodfill_highprob | awk '{print $1-0.1}'`
    echo "decreasing to floodfill_highprob=$floodfill_highprob"
    continue
endif

set test = `echo $floodfill_sigma | awk '{print ( $1 < 2 )}'`
if( $test ) then
    set floodfill_sigma = `echo $floodfill_sigma | awk '{print $1+0.1}'`
    echo "increasing to floodfill_sigma=$floodfill_sigma"
    continue
endif

set test = `echo $floodfill_Bsmooth | awk '{print ( $1 < 50 )}'`
if( $test ) then
    set floodfill_Bsmooth = `echo $floodfill_Bsmooth | awk '{print $1*2}'`
    echo "increasing to floodfill_Bsmooth=$floodfill_Bsmooth"
    continue
endif

set test = `echo $floodfill_rethresh | awk '{print ( $1 < 0.9 )}'`
if( $test ) then
    set floodfill_rethresh = `echo $floodfill_rethresh | awk '{print sqrt($1)}'`
    echo "increasing to floodfill_rethresh=$floodfill_rethresh"
    continue
endif

set floodfill_Bsmooth = `echo $floodfill_Bsmooth | awk '{print $1*2}'`
echo "incerasing smoothing B factor to $floodfill_Bsmooth"

set test = `echo $floodfill_Bsmooth | awk '{print ( $1 > 999 )}'`
if( $test ) then
  set BAD = "unable to find foothills of peak in floodfill"
  goto exit
endif

end

echo "final settings: "
cat << EOF 
floodfill_sigma=$floodfill_sigma
floodfill_Bsmooth=$floodfill_Bsmooth
floodfill_maskBsmooth=$floodfill_maskBsmooth
floodfill_rethresh=$floodfill_rethresh
floodfill_maxradius=$floodfill_maxradius
floodfill_recycles=$floodfill_recycles
floodfill_fftmargin=$floodfill_fftmargin
floodfill_highprob=$floodfill_highprob
floodfill_lowprob=$floodfill_lowprob
EOF
echo ""

mapmask mapin ${t}thresh.map mapout floodfill_mask.map << EOF >! ${t}floodmask.log
xyzlim $mapLIMITS
axis $mapAXIS
pad 0
EOF

# maybe fill in any 1-voxel holes
#map_func.com -func maxradius -param 1 floodfill_mask.map -outfile filledin.map


echo "merging floodfill_mask with significant difference mask"
mv significant_solvent_mask.map prev_significant_solvent_mask.map
echo maps mult |\
mapmask mapin1 floodfill_mask.map mapin2 prev_significant_solvent_mask.map mapout significant_solvent_mask.map >> $logfile



skip_floodfill:


set psolvent = solvent
if( ! $exclude_protein ) set psolvent = "density"
echo "fofc_signif_sol is significant ${psolvent} features only from fofc${fofc_style}"
map_func.com -func multiply fofc${fofc_style}.map significant_solvent_mask.map -outfile fofc_signif_sol.map >> $logfile
echo | mapdump mapin fofc_signif_sol.map | egrep -i "density"

echo "fofc_squish_sol is fofc with significant ${psolvent} features suppressed"
map_func.com -func negate significant_solvent_mask.map -outfile squish_solvent_mask.map >> $logfile
map_func.com -func multiply fofc${fofc_style}.map squish_solvent_mask.map -outfile fofc_squish_sol.map >> $logfile




echo "making mtz versions of these"
gemmi map2sf -v --dmin=$reso fofc_squish_sol.map gemmi.mtz dF PHIdF  >>& $logfile
echo labin file 1 all | cad hklin1 gemmi.mtz hklout fofc_squish_sol.mtz  >> $logfile
gemmi map2sf -v --dmin=$reso fofc_signif_sol.map gemmi.mtz Fpart PHIpart  >>& $logfile
echo labin file 1 all | cad hklin1 gemmi.mtz hklout fofc_signif_sol.mtz  >> $logfile
# and apply the overshoot scale
cad hklin1 fofc_signif_sol.mtz hklout sigsoldiff_overshoot.mtz << EOF >> $logfile
labin file 1 E1=Fpart E2=PHIpart
scale file 1 $overshoot_scale
labou file 1 E1=Foshoot E2=PHIoshoot
EOF

echo "adding $FC_addback back to solvent-squished version of fofc, plus $overshoot_scale overshoot, to get new FP"
set FC = `echo $FC_addback | awk -F "," '{print $1}'`
set PHIC = `echo $FC_addback | awk -F "," '{print $2}'`
cad hklin1 $mtzfile hklout cadded.mtz << EOF >> $logfile
labin file 1 E1=$FP E2=$SIGFP E3=$FC E4=$PHIC
labou file 1 E1=FP E2=SIGFP E3=FC E4=PHIC
EOF

rm -f sftooled.mtz
sftools << EOF >> $logfile
read fofc_squish_sol.mtz
read sigsoldiff_overshoot.mtz
read cadded.mtz col FP SIGFP FC PHIC
calc ( COL Fms PHI ) = ( COL FC PHIC ) ( COL dF PHIdF ) +
calc ( COL Fms PHI ) = ( COL Fms PHI ) ( COL Foshoot PHIoshoot ) +
calc F col FP = col Fms
absent col FP if col SIGFP absent
select col SIGFP = PRESENT
purge nodata yes
select all
write sftooled.mtz col FP SIGFP PHI
quit
y
EOF
echo "new FP with $FreeR_flag from $mtzfile as: $outfile"
cad hklin1 sftooled.mtz hklin2 $mtzfile hklout $outfile << EOF >> $logfile
labin file 1 E1=FP E2=SIGFP
labin file 2 E1=$FreeR_flag
EOF

echo "partial structure Fpart,PHIpart from fofc_signif_sol and $FP,$SIGFP from $mtzfile as: $refmacoutfile"
cad hklin1 $mtzfile hklin2 fofc_signif_sol.mtz hklout $refmacoutfile << EOF >> $logfile
labin file 1 E1=$FP E2=$SIGFP E3=$FreeR_flag
labin file 2 E1=Fpart E2=PHIpart
EOF

fft hklin sftooled.mtz mapout ${t}ffted.map << EOF >> $logfile
labin F1=FP PHI=PHI
$mapGRID
EOF
mapmask mapin ${t}ffted.map mapout new_FP.map << EOF >> $logfile
xyzlim $mapLIMITS
axis $mapAXIS
EOF


exit:

if("$tempfile" == "") set  tempfile = "./"
if(! $debug && ! ( "$tempfile" == "./" ) ) then
    rm -f ${tempfile}*
else
    echo "not cleaning up:"
    ls -lrt ${tempfile}*
endif

if($?BAD) then
   echo "ERROR: $BAD"
   exit 9
endif
exit

#######################################################################################
#
#  notes and example usage
#

