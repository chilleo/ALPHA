from RunDGEN import *

k = 0

#8 tax - sites ignored? (was just using this line for tsting, prolly shouldnt rerun)
#run_saved_dgen('/Users/leo/rice/res/data/dgen/tmp/figBig8/stat_8tax.txt', ['/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim1/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim2/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim3/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim4/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim5/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim6/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim7/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim8/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim9/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim10/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim11/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim12/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim13/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim14/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim15/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim16/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim17/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim18/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim19/seqfileNamed','/Users/leo/rice/res/data/dgen/simulations/12taxa/bigNetwork/firstTry/sim20/seqfileNamed'], verbose=True, plot='/Users/leo/rice/res/data/dgen/tmp/figBig8/plot_big_8', meta='Dgen')

# cichlid alignment - i put the commands for these in the java code
#run_saved_dgen('/Users/leo/rice/res/data/dgen/tmp/figC/stat_cichlid.txt', ['/Users/leo/rice/res/data/cichlid/alignment/cichlid6tax.phylip-sequential.txt'], verbose=False, plot='/Users/leo/rice/res/data/dgen/tmp/figC/plot_figC', meta='Dgen')
run_saved_dgen('/Users/leo/rice/res/data/dgen/tmp/figC/stat_cichlid.txt', ['/Users/leo/rice/res/data/cichlid/alignment/cichlid6tax.phylip-sequential.txt'], verbose=True, plot='/Users/leo/rice/res/data/dgen/tmp/figC/plot_figCVerbose', meta='Dgen')
run_saved_dgen('/Users/leo/rice/res/data/dgen/tmp/figC/stat_cichlid.txt', ['/Users/leo/rice/res/data/cichlid/alignment/cichlid6tax.phylip-sequential.txt'], window_size=10000, window_offset=500, verbose=True,  plot='/Users/leo/rice/res/data/dgen/tmp/figC/plot_figCwindows', meta='Dgen')


# Species Tree: (((P1,P2),(P3,P4)),O);
# Reticulations: [('P1', 'P3')]
test1 = generate_network_tree((0.1, 0.9), "(((P1,P2),(P3,P4)),O);", [('P1', 'P3')])
test1new = Create_Network_Helper("(((P1,P2),(P3,P4)),O);", [('P1', 'P3')],0.9)
print Create_Network_Helper("(((P1,P2),(P3,P4)),O);", [('P1', 'P3')],0.9)

#for monday, do some sims (and rerun on some old ones prolly), generate some network string here, and then run over some sims and start looking over results

run_saved_dgen("/Users/leo/rice/res/data/dgen/tmp/testStatFoil2and3iThinkO.txt",
               ['/Users/leo/rice/res/data/dgen/simulations/5taxa/noGeneFlow/dFoilStandard50kbp/sim1/seqfile','/Users/leo/rice/res/data/dgen/simulations/5taxa/noGeneFlow/dFoilStandard50kbp/sim2/seqfile'], plot='/Users/leo/rice/res/data/dgen/tmp/plot_dgen_dFoilStandard50kbp', meta='0.0')


run_saved_dgen("/Users/leo/rice/res/data/dgen/tmp/testStatFoil2and3iThinkO.txt",
               "/Users/leo/rice/res/data/dgen/simulations/5taxa/noGeneFlow/dFoilStandard50kbp/sim1/seqfile", plot='/Users/leo/rice/res/data/dgen/tmp/plot_dgen_dFoilStandard50kbp', meta='0.0')



#run_saved_dgen("/Users/leo/rice/res/data/dgen/tmp/testStatFoil2and3iThink.txt",
#               "/Users/leo/rice/res/data/dgen/tmp/M0p01")



