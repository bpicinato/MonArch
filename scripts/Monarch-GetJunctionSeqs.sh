#!/bin/bash

OUTPUT="$2"
REF_GENOME="$3"
READS="$4"
MIN_SEC_HIT_SIZE="$5"
INVERT_STRAND="$6"

while read text; do
        # Thanks LinuxHint: https://linuxhint.com/bash_split_examples/
        delimiter="franceschini-santos"
        #Concatenate the delimiter with the main string
        #Split the text based on the delimiter
        myarray=()
        while [[ $text ]]; do
                myarray+=( "${text%%"$delimiter"*}" )
                text=${text#*"$delimiter"}
        done

        # Transforming in readable variables
        CH=${myarray[0]}
        RN=${myarray[1]}
        S=${myarray[2]}
        H=${myarray[3]}
        SHC=${myarray[4]}
        OV=${myarray[5]}
        RS1=${myarray[6]}
        RE1=${myarray[7]}
        RS2=${myarray[8]}
        RE2=${myarray[9]}
        GS1=${myarray[10]}
        GE1=${myarray[11]}
        GS2=${myarray[12]}
        GE2=${myarray[13]}
        # Create tmp file with the name of the read, containig this infos
        mkdir -p ${OUTPUT}/tmp
        touch ${OUTPUT}/tmp/${RN/\//_}
        #touch $RN
        ReadSeq=$(grep -A1 "$RN" "$READS" | grep -v ">")
        # Defining vars for both cases (overlap and colinha)
        # Start coordenate of the transcript
        TRANS_START=$(echo -e "$GS1\n$GE2" | sort -n | sed '1q;d')
        # End coord. of transcript
        TRANS_END=$(echo -e "$GS1\n$GE2" | sort -n | sed '2q;d')
        # IF INVERT_STRAND==YES:::
        S_FIXD=$(echo $S | tr '+-' '-m' | tr 'm' '+')
        # Dist between half A and halfB (to check if has colinha or overlap)
        DISTANCE_BETWEEN_AnB=$(( "$RS2"-"$RE1" ))
        # Size of overlap (absolute value)
        OV_SIZE=$(echo $OV | sed 's/-//g')
        # Size of the circle
        dif=$((GE1 - GS2))
        CIRC_SIZE=$((${dif#-} + 1))
        #echo""; echo $CIRC_SIZE
        # Check if has overlap or 'colinha'
        if [[ 0 -lt ${DISTANCE_BETWEEN_AnB} ]] && [ "${RE1}" != "${RS2}" ]; then
                #echo "TEM COLINHA"; continue
                # sintaxe:
                #   ${OriginalString : FROM : LENGTH }
                READ_SEQ_A=${ReadSeq:$RS1-1:$RE1-$RS1+1}
                READ_SEQ_B=${ReadSeq:$RS2-1:$RE2-$RS2+1}
                # Size of alignment of the read
                ALI_SIZE=$((${#READ_SEQ_A} + ${#READ_SEQ_B}))
                #echo -e "READ_SEQ_A\n$READ_SEQ_A\n\nREAD_SEQ_B\n$READ_SEQ_B\n"
                COLINHA=$(echo ${ReadSeq:$RE1:$RS2-$RE1-1} | tr 'ACTG' 'actg')
                #REMAINING_TRANSCRIPT=$(echo -e "$CH\t$((GE1-1))\t$GS2\tx\t1\t$S" | \
                REMAINING_TRANSCRIPT=$(echo -e "$CH\t$TRANS_START\t$((TRANS_END-1))\tx\t1\t$S" | \
                        $CONDA_PREFIX/bin/bedtools getfasta \
                        -fi $REF_GENOME \
                        -nameOnly -bed stdin -s 2>/dev/null | \
                        grep -v ">")
                READ_SEQ_A_FXD=$(echo $READ_SEQ_A | tr ACGTacgt TGCAtgca | rev )
                READ_SEQ_B_FXD=$(echo $READ_SEQ_B | tr ACGTacgt TGCAtgca | rev )
                COLINHA_FXD=$(echo $COLINHA | tr ACGTacgt TGCAtgca | rev )
                # Only print small junction if size of the halfs are smaller
                # than 10 (ideal junction size)
                if [[ 10 -le ${#READ_SEQ_A} ]] && [[ 10 -le ${#READ_SEQ_B} ]]; then
                        READ_A_JUN=${READ_SEQ_A: -10}
                        READ_B_JUN=${READ_SEQ_B:0:10}
                else
                        READ_A_JUN=${READ_SEQ_A: -$MIN_SEC_HIT_SIZE}
                        READ_B_JUN=${READ_SEQ_B:0:$MIN_SEC_HIT_SIZE}
                fi
                READ_A_JUN_FXD=$(echo $READ_A_JUN | tr ACGTacgt TGCAtgca | rev )
                READ_B_JUN_FXD=$(echo $READ_B_JUN | tr ACGTacgt TGCAtgca | rev )

                EMPIRICAL_JUNCTION="${READ_A_JUN}|$COLINHA|${READ_B_JUN}"

                if [ "$INVERT_STRAND" == "YES" ]; then
                        # for cases in which the circle is smaller than the read::
                        if [ "$CIRC_SIZE" -lt "$ALI_SIZE" ]; then
                                READ_SEQ_B_FXD=$(echo \
                                        ${READ_SEQ_B:0:$((CIRC_SIZE-RS2-OV_SIZE+1))} \
                                        | tr ACGTacgt TGCAtgca | rev )
                                TRANSCRIPT="${COLINHA_FXD}|${READ_SEQ_A_FXD}${READ_SEQ_B_FXD}|"
                        else
                                TRANSCRIPT="${COLINHA_FXD}|${READ_SEQ_A_FXD}${REMAINING_TRANSCRIPT}${READ_SEQ_B_FXD}|"
                        fi
                        REAL_JUNCTION="${READ_B_JUN_FXD}|$COLINHA_FXD|${READ_A_JUN_FXD}"
                else
                        TRANSCRIPT="${COLINHA}|${READ_SEQ_B}${REMAINING_TRANSCRIPT}${READ_SEQ_A}|"
                        REAL_JUNCTION="${READ_A_JUN}|$COLINHA|${READ_B_JUN}"
                fi
                # Col number zero: circ size (less two, because of the || pipes)
                echo -ne "$((${#TRANSCRIPT} - 2))\t" >> ${OUTPUT}/tmp/${RN/\//_}
                # First col: first half of read
                #---------------------------------------------------------------
                echo -ne "${READ_SEQ_A}\t" >> ${OUTPUT}/tmp/${RN/\//_}
                # Second col: second half of read
                #---------------------------------------------------------------
                echo -ne "${READ_SEQ_B}\t" >> ${OUTPUT}/tmp/${RN/\//_}
                # Third col: empirical junction (Last10A..ov|colinha..First10B)
                #---------------------------------------------------------------
                echo -en "${EMPIRICAL_JUNCTION}\t" >> ${OUTPUT}/tmp/${RN/\//_}
                # Forth col: real junction (rev compl of above)
                #### IF INVERT STRAND:::::::::::::
                #---------------------------------------------------------------
                echo -en "${REAL_JUNCTION}\t" >> ${OUTPUT}/tmp/${RN/\//_}
                # Fifth col: transcript (..ov|co..ReadB..transc..ReadA)
                # Representation B'
                #---------------------------------------------------------------
                echo -en  "...${TRANSCRIPT}...\n">> ${OUTPUT}/tmp/${RN/\//_}

        else
                #echo "tem overlap"; continue
                READ_SEQ_A=${ReadSeq:$RS1-1:$RE1-$OV_SIZE-$RS1+1}
                READ_SEQ_B=${ReadSeq:$RS2-1+$OV_SIZE:$RE2-$OV_SIZE-$RS2+1}
                # Size of alignment of the read
                ALI_SIZE=$((${#READ_SEQ_A} + ${#READ_SEQ_B}))
                #echo -e "READ_SEQ_A\n$READ_SEQ_A\n\nREAD_SEQ_B\n$READ_SEQ_B\n"
                OVERLAP=$(echo ${ReadSeq:$RS2-1:$RE1-$RS2+1} | tr 'ACTG' 'actg')
                REMAINING_TRANSCRIPT=$(echo -e "$CH\t$TRANS_START\t$((TRANS_END-1))\tx\t1\t$S" |
                        $CONDA_PREFIX/bin/bedtools getfasta \
                        -fi $REF_GENOME \
                        -nameOnly -bed stdin -s 2>/dev/null | \
                        grep -v ">")
                READ_SEQ_A_FXD=$(echo $READ_SEQ_A | tr ACGTacgt TGCAtgca | rev )
                READ_SEQ_B_FXD=$(echo $READ_SEQ_B | tr ACGTacgt TGCAtgca | rev )
                OVERLAP_FXD=$(echo $OVERLAP | tr ACGTacgt TGCAtgca | rev )
                # Only print small junction if size of the halfs are smaller
                # than 10 (ideal junction size)
                if [[ 10 -le ${#READ_SEQ_A} ]] && [[ 10 -le ${#READ_SEQ_B} ]]; then
                        READ_A_JUN=${READ_SEQ_A: -$((10-OV_SIZE))}
                        READ_B_JUN=${READ_SEQ_B:0:$((10-OV_SIZE))}
                else
                        READ_A_JUN=${READ_SEQ_A: -$((MIN_SEC_HIT_SIZE-OV_SIZE))}
                        READ_B_JUN=${READ_SEQ_B:0:$((MIN_SEC_HIT_SIZE-OV_SIZE))}
                fi
                EMPIRICAL_JUNCTION="${READ_A_JUN}$OVERLAP${READ_B_JUN}"

                READ_A_JUN_FXD=$(echo $READ_A_JUN | tr ACGTacgt TGCAtgca | rev )
                READ_B_JUN_FXD=$(echo $READ_B_JUN | tr ACGTacgt TGCAtgca | rev )
                # # for cases in which the circle is smaller than the read::
                # if [ "$CIRC_SIZE" <= "$ALI_SIZE" ]; then
                #         READ_SEQ_B_FXD=$(echo ${READ_SEQ_B:0:$((CIRC_SIZE-RS2+1))} \
                #                 | tr ACGTacgt TGCAtgca | rev )
                #         TRANSCRIPT="${OVERLAP}${READ_SEQ_A_FXD}${READ_SEQ_B_FXD}"
                # fi

                if [ "$INVERT_STRAND" == "YES" ]; then
                        # for cases in which the circle is smaller than the read::
                        if [ "$CIRC_SIZE" -lt "$ALI_SIZE" ]; then
                                READ_SEQ_B_FXD=$(echo \
                                        ${READ_SEQ_B:0:$((CIRC_SIZE-RS2-OV_SIZE+1))} \
                                        | tr ACGTacgt TGCAtgca | rev )
                                TRANSCRIPT="${OVERLAP_FXD}${READ_SEQ_A_FXD}${READ_SEQ_B_FXD}"
                        else
                                TRANSCRIPT="${OVERLAP_FXD}${READ_SEQ_A_FXD}${REMAINING_TRANSCRIPT}${READ_SEQ_B_FXD}"
                        fi
                        REAL_JUNCTION="${READ_B_JUN_FXD}$OVERLAP_FXD${READ_A_JUN_FXD}"
                else
                        TRANSCRIPT="${OVERLAP}${READ_SEQ_B}${REMAINING_TRANSCRIPT}${READ_SEQ_A}"
                        REAL_JUNCTION="${READ_A_JUN}$OVERLAP${READ_B_JUN}"
                fi
                #echo -e "READ_SEQ_A\n$READ_SEQ_A$OVERLAP\n\nREAD_SEQ_B\n$OVERLAP$READ_SEQ_B\n"
                #echo -e "OVERLAP\n$OVERLAP\n"
                # Col number zero: circ size
                echo -ne "${#TRANSCRIPT}\t" >> ${OUTPUT}/tmp/${RN/\//_}
                # First col: first half of read
                #---------------------------------------------------------------
                echo -ne "${READ_SEQ_A}${OVERLAP}\t" >> ${OUTPUT}/tmp/${RN/\//_}
                # Second col: second half of read
                #---------------------------------------------------------------
                echo -ne "${OVERLAP}${READ_SEQ_B}\t" >> ${OUTPUT}/tmp/${RN/\//_}
                # Third col: empirical junction (Last10A..ov|colinha..First10B)
                #---------------------------------------------------------------
                echo -en "${EMPIRICAL_JUNCTION}\t" >> ${OUTPUT}/tmp/${RN/\//_}
                # Forth col: real junction (rev compl of above)
                #### IF INVERT STRAND:::::::::::::
                #---------------------------------------------------------------
                echo -en "${REAL_JUNCTION}\t" >> ${OUTPUT}/tmp/${RN/\//_}
                # Fifth col: transcript (..ov|co..ReadB..transc..ReadA)
                # Representation B'
                echo -en "...${TRANSCRIPT}...\n" >> ${OUTPUT}/tmp/${RN/\//_}
        fi

done <<< "$1"
