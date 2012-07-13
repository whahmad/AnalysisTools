#! /usr/bin/perl
use Cwd;
use POSIX;
use POSIX qw(strftime);

#############################################
$numArgs = $#ARGV +1;
$ARGV[$argnum];

$UserID= POSIX::cuserid();
#Default vaules
$InputDir="/net/scratch_cms/institut_3b/$UserID/Test";
$OutputDir="/net/scratch_cms/institut_3b/$UserID";
#$OutputDir="~/";
$CodeDir="../Code";
$set="TriggerStudy_";
$CMSSWRel="4_4_4";
$PileupVersion="V08-03-17";
$tag="HEAD";
$TauReco="4_4_Y_08_03_2012";
$BTag="YES";
$Cleaning ="YES";
$maxdata=100;
$maxmc=100;

if($ARGV[0] eq "--help" || $ARGV[0] eq ""){
    printf("\nThis code requires one input option. The systax is:./todo_Grid.pl [OPTION]");
    printf("\nPlease choose from the following options:\n");
    printf("\n./todo.pl --help                                   Prints this message\n");
    printf("\n./todo.pl --TauNtuple <TauNtupleDir>               Setups up TauNtuple and gives instructions for submission");
    printf("\n                                                   <TauNtupleDir> location of CMSSW for TauNtuple.");
    printf("\n                                                   Optional Commmads: ");
    printf("\n                                                     --CMSSWRel <CMSSW release #> The CMSSW release you want to use Default: $CMSSWRel");
    printf("\n                                                     --PileupVersion <Version # ie V08-03-17>");
    printf("\n                                                     --tag <tag for TauNuptle/SkimProduction/TiggerFiliter/...>");
    printf("\n                                                     --TauReco <tag for Recommended TauReco> current option: na or 4_4_Y_08_03_2012");
    printf("\n                                                     --BTag <YES/NO> Default: YES");
    printf("\n./todo.pl --Local <Input.txt>                      INTENTED FOR SMALL SCALE TESTS ONLY");  
    printf("\n                                                   Configure a directory to run locally. <InputPar.txt> name of file that");
    printf("\n                                                   contains input command template.");
    printf("\n                                                   Optional commands:  ");
    printf("\n                                                     --InputDir <InputDir>   Default value: $InputDir"); 
    printf("\n                                                     --OutputDir <OutputDir> Default value: $OutputDir");
    printf("\n                                                     --CodeDir  <CodeDir>    Default value: $CodeDir");
    printf("\n                                                     --SetName <SetName>     Default value: $set ");
    printf("\n                                                     --NMaxData <Max Number of data files per job >     Default value: $maxdata ");
    printf("\n                                                     --NMaxMC <Max Number of MC files per job >     Default value: $maxmc ");
    printf("\n./todo.pl --DCache <Input.txt> <ListofDS.txt>      INTENTED FOR REGULAR USE (DEFAULT)");
    printf("\n                                                   Configure a directory to run from. <InputPar.txt> name of file that");
    printf("\n                                                   contains input command template.");
    printf("\n                                                   <ListoDS.txt> list of DCache Dataset directories you want to run on.");
    printf("\n                                                   Optional commands:  ");
    printf("\n                                                     --OutputDir <OutputDir> ");
    printf("\n                                                     --CodeDir <CodeDir>");
    printf("\n                                                     --SetName <SetName> "); 
    printf("\n                                                     --NMaxData <Max Number of data files per job >     Default value: $maxdata ");
    printf("\n                                                     --NMaxMC <Max Number of MC files per job >     Default value: $maxmc ");
    printf("\n  ");
    exit(0);  
} 

######################################
$InputFile=$ARGV[1];
for($l=2;$l<$numArgs; $l++){
    if($ARGV[$l] eq "--InputDir"){
	$l++;
	$InputDir=$ARGV[$l];
    }
    if($ARGV[$l] eq "--OutputDir"){
	$l++;
	$OutputDir=$ARGV[$l];
    }
    if($ARGV[$l] eq "--CodeDir"){
	$l++;
	$CodeDir=$ARGV[$l];
    }
    if($ARGV[$l] eq "--SetName"){
	$l++;
	$set=$ARGV[$l];
    }
    if($ARGV[$l] eq "--CMSSWRel"){
        $l++;
        $CMSSWRel=$ARGV[$l];
    }
    if($ARGV[$l] eq "--PileupVersion"){
       $l++;
        $PileupVersion=$ARGV[$l];
    }
    if($ARGV[$l] eq "--tag"){
        $l++;
        $tag=$ARGV[$l];
    }
    if($ARGV[$l] eq "--TauReco"){
        $l++;
        $TauReco=$ARGV[$l];
    }
    if($ARGV[$l] eq "--BTag"){
	$l++;
	$BTag=$ARGV[$l];
    }
    if($ARGV[$l] eq "--Cleaning"){
	$l++;
	$Cleaning =$ARGV[$l];
    }
    if($ARGV[$l] eq  "--NMaxData" ){
	$l++;
	$maxdata=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--NMaxMC" ){
	$l++;
	$maxmc=$ARGV[$l];
    }
}

$time= strftime("%h_%d_%Y",localtime);
$temp= $set . $time;
$set=$temp;

if( $ARGV[0] eq "--TauNtuple"){
    $currentdir=getcwd;
    if($ARGV[1] ne ""){
	$basedir=$ARGV[1];
    }
    else{
	printf("\nWorkingDir for CMSSW is required. Please follow the syntax:./todo.pl --TauNtuple <TauNtupleDir> ");
	printf("\nFor more details use: ./todo --help\n"); 
	exit(0);
    }

    printf("\nWorkingDir for CMSSW: $basedir");
    printf("\nCurrentDir is: $currentdir");
    printf("\nUsing CMSSW Release: $CMSSWRel");
    printf("\nUsing PileupVersion: $PileupVersion");

    # setup CMSSW
    system(sprintf("rm Setup_TauNtuple_$CMSSWRel-$time; rm Setup_$CMSSWRel-$time "));
    system(sprintf("echo \"rm $basedir/TauNtuple_$CMSSWRel-$time -rf; mkdir $basedir/TauNtuple_$CMSSWRel-$time\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cd $basedir/TauNtuple_$CMSSWRel-$time\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    if($CMSSWRel eq "4_4_4"){
	system(sprintf("echo \"export SCRAM_ARCH=\\\"slc5_amd64_gcc434\\\"\" >> Setup_TauNtuple_$CMSSWRel-$time"));
	system(sprintf("echo \"export SCRAM_ARCH=\\\"slc5_amd64_gcc434\\\"\" >> Setup_$CMSSWRel-$time"));
    }
    if($CMSSWRel eq "5_2_0"){
        system(sprintf("echo \"export SCRAM_ARCH=\\\"slc5_amd64_gcc462\\\"\" >> Setup_TauNtuple_$CMSSWRel-$time"));
        system(sprintf("echo \"export SCRAM_ARCH=\\\"slc5_amd64_gcc462\\\"\" >> Setup_$CMSSWRel-$time"));
    }

    system(sprintf("echo \"source /net/software_cms/cmsset_default.sh\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cmsrel CMSSW_$CMSSWRel\" >> Setup_TauNtuple_$CMSSWRel-$time")); 
    system(sprintf("echo \"cd CMSSW_$CMSSWRel/src\" >> Setup_TauNtuple_$CMSSWRel-$time")); 
    system(sprintf("echo \"cmsenv\" >> Setup_TauNtuple_$CMSSWRel-$time")); 

    system(sprintf("echo \"cd $basedir/TauNtuple_$CMSSWRel-$time\" >> Setup_$CMSSWRel-$time"));
    system(sprintf("echo \"source /net/software_cms/cmsset_default.sh\" >> Setup_$CMSSWRel-$time"));
    system(sprintf("echo \"cd CMSSW_$CMSSWRel/src\" >> Setup_$CMSSWRel-$time"));
    system(sprintf("echo \"cmsenv\" >> Setup_$CMSSWRel-$time"));

    system(sprintf("echo \"addpkg RecoVertex/KinematicFit V02-00-06 \" >> Setup_TauNtuple_$CMSSWRel-$time")); 
    system(sprintf("echo \"addpkg RecoVertex/KinematicFitPrimitives  V02-00-02 \" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cvs co -d DataFormats -r V01-03-03 UserCode/RWTH3b/Tau/src/DataFormats\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cvs co -d RecoTauTag  -r V01-03-03 UserCode/RWTH3b/Tau/src/RecoTauTag\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cvs co -d CommonTools -r V01-03-03 UserCode/RWTH3b/Tau/src/CommonTools\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"addpkg PhysicsTools/Utilities $PileupVersion \" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"addpkg RecoLuminosity/LumiDB V03-03-16\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cvs co -d TauDataFormat  -r $tag UserCode/RWTH3b/Tau/FlatNtuple/TauDataFormat/\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cvs co -d SkimmingTools  -r $tag UserCode/RWTH3b/Tau/FlatNtuple/SkimmingTools/\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cvs co -d SkimProduction -r $tag UserCode/RWTH3b/Tau/FlatNtuple/SkimProduction/\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cvs co -d TriggerFilter  -r $tag UserCode/RWTH3b/Tau/FlatNtuple/TriggerFilter\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cp $currentdir/subs SkimProduction/CRAB/\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    if($BTag eq "YES"){
	system(sprintf("echo \"cvs co -r V00-04-11 RecoBTag/PerformanceDB\" >> Setup_TauNtuple_$CMSSWRel-$time"));
	system(sprintf("echo \"cvs co -d TopAnalysis/TopUtils UserCode/Bromo/TopAnalysis/TopUtils\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    }
    if($TauReco eq "4_4_Y_08_03_2012"){
	system(sprintf("echo \"cvs co -r V01-02-01 RecoTauTag/TauTagTools \" >> Setup_TauNtuple_$CMSSWRel-$time"));
	system(sprintf("echo \"cvs co -r V01-02-16 RecoTauTag/RecoTau \" >> Setup_TauNtuple_$CMSSWRel-$time"));
	system(sprintf("echo \"cvs co -r V01-02-12 RecoTauTag/Configuration \" >> Setup_TauNtuple_$CMSSWRel-$time"));
	system(sprintf("echo \"addpkg  PhysicsTools/PatAlgos \" >> Setup_TauNtuple_$CMSSWRel-$time"));
	system(sprintf("echo \"cvs up -r 1.47 PhysicsTools/PatAlgos/python/tools/tauTools.py \" >> Setup_TauNtuple_$CMSSWRel-$time"));
    }
    # add private hacks
    if($CMSSWRel eq "4_2_8"){    
	system(sprintf("echo \"$currentdir/subs \\\"AlgebraicMatrix33\\\" \\\"AlgebraicSymMatrix\\\" RecoVertex/KinematicFitPrimitives/src/KinematicVertex.cc\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    }
    if($PileupVersion eq "V08-03-14"){
	system(sprintf("echo \"$currentdir/subs \\\"PUInputFile_,PUInputFile_, PUInputHistoMC_, PUInputHistoData_,PUOutputFile_\\\" \\\"PUInputFile_,PUInputFile_, PUInputHistoMC_, PUInputHistoData_\\\" TauDataFormat/TauNtuple/src/TauNtuple.cc\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    }
    if($Cleaning eq "YES"){
	system(sprintf("echo \"cvs co -r lhx_12JAN2012_v1 DataFormats/METReco\" >> Setup_TauNtuple_$CMSSWRel-$time"));
	system(sprintf("echo \"cvs co -r lhx_08FEB2012_v2 RecoMET/METFilters\" >> Setup_TauNtuple_$CMSSWRel-$time"));
	system(sprintf("echo \"cvs co -r V00-00-08 RecoMET/METAnalyzers\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    }

    # Setup CRAB
    system(sprintf("echo \"export VO_CMS_SW_DIR=\\\"/net/software_cms\\\"\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"source /net/software_cms/cmsset_default.sh\" >> Setup_TauNtuple_$CMSSWRel-$time "));
    system(sprintf("echo \"source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh\" >> Setup_TauNtuple_$CMSSWRel-$time"));

    system(sprintf("echo \"export VO_CMS_SW_DIR=\\\"/net/software_cms\\\"\" >> Setup_$CMSSWRel-$time"));
    system(sprintf("echo \"source /net/software_cms/cmsset_default.sh\" >> Setup_$CMSSWRel-$time "));
    system(sprintf("echo \"source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh\" >> Setup_$CMSSWRel-$time"));

   #build
    system(sprintf("echo \"scram build -c\" >> Setup_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"scram b\" >> Setup_TauNtuple_$CMSSWRel-$time"));

    # print Instructions
    printf("\n\nInstructions\n");
    printf("Please setup CVSROOT, klog and kinit before proceeding...");
    printf("\nexport CVSROOT=\":gserver:cmscvs.cern.ch:/local/reps/CMSSW\"");
    printf("\nkinit <UserID>@CERN.CH"); 
    printf("\nklog <UserID>@CERN.CH"); 
    printf("\nsource Setup_TauNtuple_$CMSSWRel-$time"); 
    printf("\nYou can now go to the Production Directory and submit jobs.");
    printf("\ncd  $basedir/SkimProduction/CRAB/");
    printf("\n1) Read the instruction for todo.pl: ./todo.pl --help");
    printf("\n2) Setup the jobs: ./todo.pl --Submit HTL_Tau_Ntuple_cfg.py Input.dat ");
    printf("\n3) Setup the jobs and submit 1 Test job: ./todo.pl --Submit HTL_Tau_Ntuple_cfg.py Input.dat --n 1");
    printf("\n4) Setup the jobs and start production:./todo.pl --Submit HTL_Tau_Ntuple_cfg.py Input.dat --n -1\n");
    printf("\nIf you want to use this release in a new window: source Setup_$CMSSWRel-$time \n");
}


if( $ARGV[0] eq "--Local" ){
    # Print out input parameters
    printf("Active directory will be: $OutputDir/workdir$set \n");
    printf("Input directory will be:  $InputDir \n");
    printf("Code Repository is:       $CodeDir \n");
    # Clean Directory in case workdir$set exists
    printf("Cleaning Directories \n");
    system(sprintf("cd $OutputDir"));
    system(sprintf("rm -rf $OutputDir/workdir$set \n"));
    system(sprintf("mkdir $OutputDir/workdir$set "));
    printf("Cleaning complete \n");
    
    # create directory stucture
    system(sprintf("mkdir  $OutputDir/workdir$set/Code "));
    system(sprintf("mkdir  $OutputDir/workdir$set/Code/i386_linux "));
    system(sprintf("cp -r $CodeDir/* $OutputDir/workdir$set/Code/ "));
    system(sprintf("mkdir $OutputDir/workdir$set/EPS "));
    system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/InputData "));
    
    # generate compile script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"source config\"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake cleanall \" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake all \" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/compile")) ;
 
    # Generate Combine script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Combine")) ;
    system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Combine")); 
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Combine")) ; 
    system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe \" >> $OutputDir/workdir$set/Combine")) ;

    # Generate Combine Input
    system(sprintf("cp $InputFile $OutputDir/workdir$set/Input.txt "));
    system(sprintf("echo \"Mode: RECONSTRUCT\" >> $OutputDir/workdir$set/Input.txt"));
    system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Input.txt"));

    # Setup Condor Combine scripts
    system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"executable   = Combine      \"  >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"output       = Combine-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"error        = Combine-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"log          = Combine-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Condor_Combine")); 			
    system(sprintf("echo \"queue = 1 \" >> $OutputDir/workdir$set/Condor_Combine"));

    # Start Submit script
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Submit")) ; 
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"rm Set*/*.o; rm Set*/*.e; rm Set*/*.log; \" >> $OutputDir/workdir$set/Submit")) ;

    opendir(DIR,"$InputDir/");
    @dirs = grep {( /user/ || /data/ || /mc/)} readdir(DIR);
    closedir DIR;
    
    $B=0;
    $max=1;
    for($l=0;$l<2; $l++){
	printf("\n\nStarting Loop $l \n");
	$A=$maxdata+$maxmc+10;
	foreach $subdir (@dirs){
	    printf("$subdir");
	    if(($l==0 && ($subdir =~ m/data/)) || ($l==1 && !($subdir =~ m/data/))){
		if($l==0){
		    $max=$maxdata;
		}
		else{
		    $max=$maxmc;
		}
		printf("Accessing Directory  $subdir \n");
		opendir(SUBDIR,"$InputDir/$subdir/");
		@files = grep { /root/ } readdir(SUBDIR);

		$nfiles = @files;
		$idx=0;
		foreach $file (@files){
		    $idx++;
		    printf("$file Set= $B  Index=  $A   Max.= $max N Files = $nfiles Current File = $idx \n");
		    if($A > $max){
			$A=1;
			$B++;
			# Add Set information to Combining scripts and Input.txt
			system(sprintf("echo \"File: $OutputDir/workdir$set/Set_$B/  \" >>  $OutputDir/workdir$set/Input.txt ")) ;
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B \" >> $OutputDir/workdir$set/Submit")) ;
			system(sprintf("echo \"condor_submit  Condor_Set_$B  \" >> $OutputDir/workdir$set/Submit")) ;

			# Create and configure Set_$B dir
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B ")) ;
			system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/Set_$B/InputData "));
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B/EPS ")) ;

                        # Setup Set_$B.sh
			system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("echo \"echo 'Starting Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("chmod +x $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe  | tee Set_$B.output \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"echo 'Completed Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));

			# Setup Input.txt
			system(sprintf("cp   $InputFile $OutputDir/workdir$set/Set_$B/Input.txt ")); 
			system(sprintf("echo \"Mode: ANALYSIS\" >> $OutputDir/workdir$set/Set_$B/Input.txt")); 
			system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Set_$B/Input.txt"));

			# Setup Condor scripts
			system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"executable   = Set_$B.sh      \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"output       = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"error        = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"log          = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 			
			system(sprintf("echo \"queue = 1 \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
		    }
		    system(sprintf("echo \"File: $InputDir/$subdir/$file  \" >> $OutputDir/workdir$set/Set_$B/Input.txt")) ;
		    $A++;
		}
	    }
	}
    }

    # Finish Submit script
    system(sprintf("echo \"cd  $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"chmod +x  Auto_Combine_submit  \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"nohup  ./Auto_Combine_submit & \" >> $OutputDir/workdir$set/Submit"));
    GenerateAutoCombineScript("Auto_Combine_submit");
    system(sprintf("mv Auto_Combine_submit $OutputDir/workdir$set/"));

    # print Instructions
    printf("\n\nInstructions");
    printf("\nTo go to the Test workdir: cd  $OutputDir/workdir$set ");
    printf("\nTo compile the code in the workdir: source compile ");
    printf("\nTo submit jobs to the batch queue: source Submit ");
    printf("\nTo combine jobs submitted to the batch queue: source Combine \n");
    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile; cd $OutputDir/workdir$set/Set_1; source Set_1.sh | tee log; cd ..\n");
} 



if( $ARGV[0] eq "--DCache" ){

    $TempDataSetFile=$ARGV[2];
    # Print out input parameters
    printf("Active directory will be: $OutputDir/workdir$set \n");
    printf("Code Repository is:       $CodeDir \n");
    printf("List of dcache dir:       $TempDataSetFile \n");

    # Open ListofFile.txt
    @DataSets;
    open(DAT, $TempDataSetFile) || die("Could not open file $TempDataSetFile!");
    while ($item = <DAT>) {
	chomp($item);
	print("File: $item \n");
	push(@DataSets,$item);
    }
    close(DAT);

    # Clean Directory in case workdir$set exists
    printf("Cleaning Directories \n");
    system(sprintf("cd $OutputDir"));
    system(sprintf("rm -rf $OutputDir/workdir$set \n"));
    system(sprintf("mkdir $OutputDir/workdir$set "));
    printf("Cleaning complete \n");
    
    # creat directory stucture
    system(sprintf("mkdir  $OutputDir/workdir$set/Code "));
    system(sprintf("mkdir  $OutputDir/workdir$set/Code/i386_linux "));
    system(sprintf("cp -r $CodeDir/* $OutputDir/workdir$set/Code/ "));
    system(sprintf("mkdir $OutputDir/workdir$set/EPS "));
    system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/InputData "));
    
    # generate compile script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"source config\"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake cleanall \" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake all \" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/compile")) ;
 
    # Generate Combine script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Combine")) ;
    system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Combine"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Combine")) ; 
    system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe \" >> $OutputDir/workdir$set/Combine")) ;

    # Generate Combine Input
    system(sprintf("cp $InputFile $OutputDir/workdir$set/Input.txt "));
    system(sprintf("echo \"Mode: RECONSTRUCT\" >> $OutputDir/workdir$set/Input.txt"));
    system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Input.txt"));

    # Setup Condor Combine scripts
    system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"executable   = Combine      \"  >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"output       = Combine-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"error        = Combine-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"log          = Combine-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Condor_Combine")); 			
    system(sprintf("echo \"queue = 1 \" >> $OutputDir/workdir$set/Condor_Combine"));

    # Start Submit script
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Submit")) ; 
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"rm Set*/*.o; rm Set*/*.e; rm Set*/*.log; \" >> $OutputDir/workdir$set/Submit")) ;

    for($l=0;$l<2; $l++){
	$B=0;
	$max=1;
	foreach $DS (@DataSets){
	    if(($l==0 && ($DS =~ m/data/)) || ($l==1 && !($DS =~ m/data/))){
		if($l==0){
		    $max=$maxdata;
		}
		else{
		    $max=$maxmc;
		}
		printf("\n\nStarting Loop $l \n");
		$A=$maxdata+$maxmc+10;

		# find the root files for the current DataSet (DS)
		printf("Accessing Directory  $DS \n");
		system(sprintf("touch junk"));
		system(sprintf("uberftp grid-ftp.physik.rwth-aachen.de \"cd $DS; ls */  \" | grep root >& junk "));
		system(sprintf("cat junk | awk '{print \$9}' >& junk1")); 

		# Get list of files in dcache dir
		@files=();
		open(DAT, "junk1");
		while ($item = <DAT>) {
		    chomp($item);
		    push(@files,$item);
		}
		close(DAT);
		system(sprintf("rm junk"));
		system(sprintf("rm junk1"));
		$nfiles = @files;
		$idx=0;
		foreach $file (@files){
		    $idx++;
		    printf("$DS/$file Set= $B  Index=  $A   Max.= $max N Files = $nfiles Current File = $idx \n");
		    if($A > $max ){
			$A=1;
			$B++;

			# Add Set information to Combining scripts and Input.txt
			system(sprintf("echo \"File: $OutputDir/workdir$set/Set_$B/ \" >>  $OutputDir/workdir$set/Input.txt ")) ;
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B \" >> $OutputDir/workdir$set/Submit")) ;
			system(sprintf("echo \"condor_submit  Condor_Set_$B  \" >> $OutputDir/workdir$set/Submit")) ;


			# Create and configure Set_$B dir
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B ")) ;
			system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/Set_$B/InputData "));
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B/EPS ")) ;
			
			# Setup Set_$B.sh
			system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("echo \"echo 'Starting Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"source $OutputDir/workdir$set/Set_$B/Set_$B-get.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ; 
			system(sprintf("chmod +x $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"mkdir  /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cp -r *   /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe | tee Set_$B.output \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cp -r *  $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"source $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"rm -r  /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));			
			system(sprintf("echo \"echo 'Completed Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));

                        # Setup Set_$B_get.sh and Set_$B_clean.sh
			system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			system(sprintf("echo \"mkdir /user/scratch/$UserID \" >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			system(sprintf("echo \"cd /user/scratch/$UserID \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh")); 
			system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));
			system(sprintf("echo \"cd /user/scratch/$UserID \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh")); 

			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			# Setup Input.txt
			system(sprintf("cp   $InputFile $OutputDir/workdir$set/Set_$B/Input.txt ")); 
			system(sprintf("echo \"Mode: ANALYSIS\" >> $OutputDir/workdir$set/Set_$B/Input.txt")); 
			system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Set_$B/Input.txt"));
			
			# Setup Condor scripts
			system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"executable   = Set_$B.sh      \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"output       = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"error        = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"log          = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 			
			system(sprintf("echo \"queue = 1 \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
		    }
		    ($a,$b,$c)=split('/',$file);
		    $myfile=$file;
		    if($a =~ m/root/){
			$myfile=$a;
		    }
		    if($b =~ m/root/){
                        $myfile=$b;
                    }
		    if($c =~ m/root/){
                        $myfile=$c;
                    }
		    #printf("$myfile");
		    system(sprintf("echo \"/opt/d-cache/dcap/bin/dccp  dcap://grid-dcap.physik.rwth-aachen.de/$DS/$file . \"  >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
		    system(sprintf("echo \"File:  /user/scratch/$UserID/$myfile \"     >> $OutputDir/workdir$set/Set_$B/Input.txt")) ;
		    system(sprintf("echo \"rm -rf /user/scratch/$UserID/$myfile  \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));
		    $A++;
		}
	    }
	}
    }

    # Finish Submit script
    system(sprintf("echo \"cd  $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"chmod +x  Auto_Combine_submit  \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"nohup  ./Auto_Combine_submit & \" >> $OutputDir/workdir$set/Submit"));
    GenerateAutoCombineScript("Auto_Combine_submit");
    system(sprintf("mv Auto_Combine_submit $OutputDir/workdir$set/"));

    # print Instructions
    printf("\n\nInstructions");
    printf("\nPlease make sure you have run:");
    printf("\nvoms-proxy-init");
    printf("\ngrid-proxy-init");
    printf("\nkinit $UserID@CERN.CH");
    printf("\nNow you can run the analysis using dcache.");
    printf("\nTo go to the Test workdir: cd  $OutputDir/workdir$set ");
    printf("\nTo compile the code in the workdir: source compile ");
    printf("\nTo submit jobs to the batch queue: source Submit ");
    printf("\nTo combine jobs submitted to the batch queue: source Combine \n");
    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile; cd $OutputDir/workdir$set/Set_1; source Set_1 | tee log; cd ..\n");
} 



sub GenerateAutoCombineScript{
    $file=$_[0];
    system(sprintf("echo \"#! /bin/bash\" >> $file"));
    system(sprintf("echo \"for i in {1..2881} \" >>  $file"));
    system(sprintf("echo \"do \" >>  $file"));
    system(sprintf("echo \"  nDir=\\\$(ls | grep Set_ | wc -l) \" >>  $file"));
    system(sprintf("echo \"  nFin=\\\$(grep \\\"Completed Job\\\" Set_*/Set_*.o | wc -l) \" >>  $file"));
    system(sprintf("echo \"  if [ \\\$nFin != \\\$nDir ]; then\" >>  $file"));
    system(sprintf("echo \"    sleep 300\" >>  $file"));
    system(sprintf("echo \"  else\" >>  $file"));
    system(sprintf("echo \"    sleep 180\" >>  $file"));
    system(sprintf("echo \"    nSuc=\\\$(grep \\\"Program is Finished\\\" Set_*/Set_*.output  | wc -l )  \" >>  $file"));
    system(sprintf("echo \"    if [ \\\$nFin == \\\$nSuc ]; then\" >>  $file"));
    system(sprintf("echo \"      echo \\\"All Jobs Succeeded: Submitting Condor_Combine \\\"\" >>  $file"));
    system(sprintf("echo \"      condor_submit Condor_Combine\" >>  $file"));
    system(sprintf("echo \"    else \" >>  $file"));
    system(sprintf("echo \"      echo \\\"Failed Jobs: Aborting combine submission \\\"\" >>  $file"));
    system(sprintf("echo \"    fi\" >>  $file"));
    system(sprintf("echo \"    break\" >>  $file"));
    system(sprintf("echo \"  fi\" >>  $file"));
    system(sprintf("echo \"  if [ \\\$i == 2880 ]; then\" >>  $file"));
    system(sprintf("echo \"    echo \\\"Out of time (48hr): Submit Condor_Combine manually when jobs complete \\\" \" >>  $file"));
    system(sprintf("echo \"    break\" >>  $file")); 
    system(sprintf("echo \"  fi\" >>  $file")); 
    system(sprintf("echo \"done\" >>  $file")); 
}

