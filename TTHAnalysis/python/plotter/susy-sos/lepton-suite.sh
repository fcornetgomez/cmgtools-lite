################################
#  use mcEfficiencies.py to make plots of efficiencies and fake rates
################################

T="/data/peruzzi/soft_lep_training"
BCORE=" --s2v --tree treeProducerSusyMultilepton susy-sos/lepton-mca.txt object-studies/lepton-perlep.txt -P ${T} "
PBASE="~/www/lepton_suite/plots/74X/lepMVA/soft/v1.0"
BG=" -j 8 "
DEN=" LepGood_miniRelIso<0.4 && LepGood_sip3d < 8 && abs(LepGood_eta)<1.5"
if [[ "$1" == "-b" ]]; then BG=" & "; shift; fi

what=$1;
echo "mkdir -p $PBASE/$what 2> /dev/null && cp /afs/cern.ch/user/g/gpetrucc/php/index.php $PBASE/$what/"
case $what in
rocs)
    ROCS0="rocCurves.py -p TT_red,TT_true ${BCORE} susy-sos/lepton-sels.txt susy-sos/lepton-xvars.txt --splitSig 1 ";
    ROCS="${ROCS0}  --logx --grid --xrange 0.007 1.0 --yrange 0.80 1 --max-entries 5000000 "
    echo "( python $ROCS  -o $PBASE/$what/mu_pt_25_inf.root  -R pt20 pt '${DEN} && LepGood_pt > 25'  --sP 'mvaTTH.*,sosID' ${BG} )"
    ROCS="${ROCS0}  --logx --grid --xrange 0.007 1.0 --yrange 0.50 1 --max-entries 5000000 "
    echo "( python $ROCS  -o $PBASE/$what/mu_pt_10_25.root   -R pt20 pt '${DEN} && LepGood_pt < 25 && LepGood_pt > 10'  --sP 'mvaTTH.*,sosID' ${BG} )"
    ROCS="${ROCS0}  --logx --grid --xrange 0.0007 1.0 --yrange 0.00 1 --max-entries 5000000 "
    echo "( python $ROCS  -o $PBASE/$what/mu_pt_5_10.root    -R pt20 pt '${DEN} && LepGood_pt < 10 && LepGood_pt > 5'  --sP 'mvaTTH.*,sosID' ${BG} )"

    ROCS="${ROCS0}  --logx --grid --xrange 0.007 1.0 --yrange 0.80 1 --max-entries 5000000 "
    echo "( python $ROCS  -o $PBASE/$what/el_pt_25_inf.root  -R pt20 pt '${DEN} && LepGood_pt > 25' -I mu  --sP 'mvaTTH.*,sosID' ${BG} )"
    ROCS="${ROCS0}  --logx --grid --xrange 0.007 1.0 --yrange 0.50 1 --max-entries 5000000 "
    echo "( python $ROCS  -o $PBASE/$what/el_pt_10_25.root   -R pt20 pt '${DEN} && LepGood_pt < 25 && LepGood_pt > 10' -I mu  --sP 'mvaTTH.*,sosID' ${BG} )"
    ROCS="${ROCS0}  --logx --grid --xrange 0.0007 1.0 --yrange 0.00 1 --max-entries 500000 "
    echo "( python $ROCS  -o $PBASE/$what/el_pt_7_10.root    -R pt20 pt 'LepGood_pt < 10 && LepGood_pt > 7'  -I mu --sP 'mvaTTH.*,sosID' ${BG} )"
    ;;
effs)
    EFFS="mcEfficiencies.py -p TTH_true,TTH_true_tau,TTH_true_top ${BCORE} susy-sos/lepton-sels.txt susy-sos/lepton-xvars.txt --groupBy process,cut --legend=BR";
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/mu_eta_00_12_comp.root  -R pt20 eta '${DEN} && abs(LepGood_eta)<1.2' --sP pt_fine --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/mu_eta_12_24_comp.root  -R pt20 eta '${DEN} && abs(LepGood_eta)>1.2' --sP pt_fine --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/el_EB_comp.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)<1.479' --sP pt_fine --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/el_EE_comp.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)>1.479' --sP pt_fine --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0.7 1.  -o $PBASE/$what/mu_pt25_comp.root  -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fine --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0.7 1.  -o $PBASE/$what/el_pt25_comp.root -I mu -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fine --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0.0 1.  -o $PBASE/$what/mu_pt10_25_comp.root  -R pt20 eta '${DEN} && LepGood_pt < 25 && LepGood_pt > 10' --sP eta_fine --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0.0 1.  -o $PBASE/$what/el_pt10_25_comp.root -I mu -R pt20 eta '${DEN} && LepGood_pt < 25 && LepGood_pt > 10' --sP eta_fine --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/mu_eta_00_12_hist_T.root  -R pt20 eta '${DEN} && abs(LepGood_eta)<1.2' --sP pt_fine --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/mu_eta_12_24_hist_T.root  -R pt20 eta '${DEN} && abs(LepGood_eta)>1.2' --sP pt_fine --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/el_EB_hist_T.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)<1.479' --sP pt_fine --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/el_EE_hist_T.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)>1.479' --sP pt_fine --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0.7 1.  -o $PBASE/$what/mu_pt25_hist_T.root  -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fine --sP mvaTT.*_T ${BG} )"
    echo "( python $EFFS --yrange 0.7 1.  -o $PBASE/$what/el_pt25_hist_T.root -I mu -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fine --sP mvaTT.*_T ${BG} )"
    echo "( python $EFFS --yrange 0.0 1.  -o $PBASE/$what/mu_pt10_25_hist_T.root  -R pt20 eta '${DEN} && LepGood_pt < 25 && LepGood_pt > 10' --sP eta_fine --sP mvaTT.*_T ${BG} )"
    echo "( python $EFFS --yrange 0.0 1.  -o $PBASE/$what/el_pt10_25_hist_T.root -I mu -R pt20 eta '${DEN} && LepGood_pt < 25 && LepGood_pt > 10' --sP eta_fine --sP mvaTT.*_T ${BG} )"
    ;;
fakes)
    EFFS="mcEfficiencies.py -p TT_red,TT_bjets ${BCORE} susy-sos/lepton-sels.txt susy-sos/lepton-xvars.txt --groupBy process,cut --legend=TL";
    echo "( python $EFFS --yrange 0 0.12  -o $PBASE/$what/mu_eta_00_12_comp.root  -R pt20 eta '${DEN} && abs(LepGood_eta)<1.2' --sP pt_fcoarse --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 0.12  -o $PBASE/$what/mu_eta_12_24_comp.root  -R pt20 eta '${DEN} && abs(LepGood_eta)>1.2' --sP pt_fcoarse --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 0.12  -o $PBASE/$what/mu_pt25_comp.root  -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fcoarse --sP mvaTTH_[TM],multiIso ${BG} )"
    EFFS="mcEfficiencies.py -p TT_fake ${BCORE} susy-sos/lepton-sels.txt susy-sos/lepton-xvars.txt --groupBy process,cut --legend=TL";
    echo "( python $EFFS --yrange 0 0.35  -o $PBASE/$what/mu_eta_00_12_comp.root  -R pt20 eta '${DEN} && abs(LepGood_eta)<1.2' --sP pt_fvcoarse --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 0.35  -o $PBASE/$what/mu_eta_12_24_comp.root  -R pt20 eta '${DEN} && abs(LepGood_eta)>1.2' --sP pt_fvcoarse --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 0.35  -o $PBASE/$what/mu_pt25_comp.root  -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fvcoarse --sP mvaTTH_[TM],multiIso ${BG} )"
    EFFS="mcEfficiencies.py -p TT_red,TT_bjets,TT_fake ${BCORE} susy-sos/lepton-sels.txt susy-sos/lepton-xvars.txt --groupBy process,cut --legend=TL";
    echo "( python $EFFS --yrange 0 0.1  -o $PBASE/$what/el_EB_comp.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)<1.479' --sP pt_fcoarse --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 0.1  -o $PBASE/$what/el_EE_comp.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)>1.479' --sP pt_fcoarse --sP mvaTTH_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 0.1  -o $PBASE/$what/el_pt25_comp.root -I mu -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fcoarse --sP mvaTTH_[TM],multiIso ${BG} )"

    EFFS="mcEfficiencies.py -p TT_red,TT_bjets ${BCORE} susy-sos/lepton-sels.txt susy-sos/lepton-xvars.txt --groupBy process,cut --legend=TL";
    echo "( python $EFFS --yrange 0 0.12  -o $PBASE/$what/mu_eta_00_12_hist_T.root  -R pt20 eta '${DEN} && abs(LepGood_eta)<1.2' --sP pt_fcoarse --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 0.12  -o $PBASE/$what/mu_eta_12_24_hist_T.root  -R pt20 eta '${DEN} && abs(LepGood_eta)>1.2' --sP pt_fcoarse --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 0.12  -o $PBASE/$what/mu_pt25_hist_T.root  -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fcoarse --sP mvaTTH.*_T ${BG} )"
    EFFS="mcEfficiencies.py -p TT_fake ${BCORE} susy-sos/lepton-sels.txt susy-sos/lepton-xvars.txt --groupBy process,cut --legend=TL";
    echo "( python $EFFS --yrange 0 0.35  -o $PBASE/$what/mu_eta_00_12_hist_T.root  -R pt20 eta '${DEN} && abs(LepGood_eta)<1.2' --sP pt_fvcoarse --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 0.35  -o $PBASE/$what/mu_eta_12_24_hist_T.root  -R pt20 eta '${DEN} && abs(LepGood_eta)>1.2' --sP pt_fvcoarse --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 0.35  -o $PBASE/$what/mu_pt25_hist_T.root  -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fvcoarse --sP mvaTTH.*_T ${BG} )"
    EFFS="mcEfficiencies.py -p TT_red,TT_bjets,TT_fake ${BCORE} susy-sos/lepton-sels.txt susy-sos/lepton-xvars.txt --groupBy process,cut --legend=TL";
    echo "( python $EFFS --yrange 0 0.1  -o $PBASE/$what/el_EB_hist_T.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)<1.479' --sP pt_fcoarse --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 0.1  -o $PBASE/$what/el_EE_hist_T.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)>1.479' --sP pt_fcoarse --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 0.1  -o $PBASE/$what/el_pt25_hist_T.root -I mu -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fcoarse --sP mvaTTH.*_T ${BG} )"

    ;;
esac

