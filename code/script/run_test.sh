EXEC_PATH=/Users/forero/Dropbox/SparseFilaments/code/ngl-beta/binsrc/
DATA_IN=../../data/input/
DATA_OUT=../../data/output/
POS_FILE=Pos_MD_full_halos_900-1000kms_z0.dat
METHOD=BSkeleton
$EXEC_PATH./getNeighborGraph -i $DATA_IN$POS_FILE -d 3 -b 1.5 -m $METHOD > $DATA_OUT$METHOD\_b1.5_$POS_FILE
