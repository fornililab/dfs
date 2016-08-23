######################################################################
# T01: dfs
cd T01
bash T01.sh
cd ..
cksum T01/fixed_*_scores_T01.dat > check/T01.out
diff check/T01.out check/T01.chk > check/T01.diff
if [ -s check/T01.diff ]
then
        echo "TEST T01: FAILED"
else
        echo "TEST T01: PASSED"
fi
