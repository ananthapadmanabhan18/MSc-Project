#!/bin/bash


####################
#    Parameters    #
####################
run="294925"
mode="sim"
detector="FOCAL"
uid="1"
simulation="NoDigitization"
geometry_location="GeometryFiles/geometry.txt"
nevents="10"
generator="Upgrade:FOCAL_Generators:boxMomUniform"
pdg="111"
etamin="3.4"
etamax="5.3"
ptmin="0.0"
ptmax="1500.0"

num=$1
for ((i=1; i<=num; i++))
do
  echo "#############################"
  echo "#  Running simulation: $i   #"
  echo "#############################"





  ####################
  #    Target_DIR    #
  ####################
  SOURCE_DIR="."
  TARGET_DIR="../Simulation_outputs/pdg"$pdg"_nevents$nevents/$i"
  mkdir -p $TARGET_DIR

  ##########################
  #    Print_parameters    #
  ##########################
  echo "Running simulation with the following parameters:"
  echo "Source directory      :       $SOURCE_DIR"
  echo "Target directory      :       $TARGET_DIR"
  echo "UID                   :       $uid"
  echo "mode                  :       $mode"
  echo "simulation            :       $simulation"
  echo "Run number            :       $run"
  echo "Mode                  :       $mode"
  echo "Detector              :       $detector"
  echo "Generator             :       $generator"
  echo "PDG ID                :       $pdg"
  echo "Number of events      :       $nevents"
  echo "Eta range             :       $etamin - $etamax"
  echo "Pt range              :       $ptmin - $ptmax"


  ##########################
  #     Run Simulation     #
  ##########################
  source /home/apmnair18/alice/sw/ubuntu2204_x86-64/AliDPG/master-local1/bin/aliroot_dpgsim.sh --run $run --mode $mode --detector $detector --uid $uid --nevents $nevents --generator $generator --simulation $simulation --focalGeometryFile $geometry_location --pdg $pdg --etamin $etamin --etamax $etamax --ptmin $ptmin --ptmax $ptmax --jobs 20 --debug

  ###############################################
  # Move the GRP folder to the target directory #
  ###############################################
  if [ -d "$SOURCE_DIR/GRP" ]; then
    mv "$SOURCE_DIR/GRP" "$TARGET_DIR/"
  fi


  ######################################################
  # Move all output ROOT files to the target directory #
  ######################################################
  for file in "$SOURCE_DIR"/*.root; do
    if [[ "$file" != "$SOURCE_DIR/OCDBsim.root" && "$file" != "$SOURCE_DIR/OCDBrec.root" ]]; then
      mv "$file" "$TARGET_DIR/"
    fi
  done
  mv "$SOURCE_DIR/grpdump.sh" "$TARGET_DIR/"
  mv "$SOURCE_DIR/gphysi.dat" "$TARGET_DIR/"
  mv "$SOURCE_DIR/sim.log" "$TARGET_DIR/"
  mv "$SOURCE_DIR/simwatch.log" "$TARGET_DIR/"
  mv "$SOURCE_DIR/MCStepLoggerVolMap.dat" "$TARGET_DIR/"
done






