#include "FFTtools.h"
#include "Analyzer.h"
#include "FilteredAnitaEvent.h"
#include "AnitaDataset.h" 
#include "RawAnitaHeader.h"
#include "AnalysisConfig.h"
#include "UCFilters.h"
#include "BasicFilters.h"
#include "AnitaTemplates.h"
//#include "Util.h"
#include "FilterStrategy.h"
#include "TFile.h"
#include "TTree.h"

typedef struct {

	Float_t theta0;
	Float_t theta1;
	Float_t theta2;
	
	Float_t phi0;
	Float_t phi1;
	Float_t phi2;

	Float_t mapPeak0;
	Float_t mapPeak1;
	Float_t mapPeak2;

	Float_t mapRMS0;
	Float_t mapRMS1;
	Float_t mapRMS2;

	Float_t coherent_hpeak0;
	Float_t coherent_hpeak1;
	Float_t coherent_hpeak2;
	Float_t coherent_snr0;	
	Float_t coherent_snr1;	
	Float_t coherent_snr2;	
	Float_t coherent_filtered_hpeak0;
	Float_t coherent_filtered_hpeak1;
	Float_t coherent_filtered_hpeak2;
	Float_t coherent_filtered_snr0;	
	Float_t coherent_filtered_snr1;	
	Float_t coherent_filtered_snr2;	

	Float_t deconv_hpeak0;
	Float_t deconv_hpeak1;
	Float_t deconv_hpeak2;
	Float_t deconv_snr0;	
	Float_t deconv_snr1;	
	Float_t deconv_snr2;	
	Float_t deconv_filtered_hpeak0;
	Float_t deconv_filtered_hpeak1;
	Float_t deconv_filtered_hpeak2;
	Float_t deconv_filtered_snr0;	
	Float_t deconv_filtered_snr1;	
	Float_t deconv_filtered_snr2;	

} hPol,vPol;


typedef struct {

	Int_t eventNumber;
	Int_t run;
//	Int_t isWais;
	Int_t isRF;
	Int_t isPayloadBlast;
	Int_t isVarner;

} eventInfo;

typedef struct {

	Float_t maxBottomToTopRatioH;
	Int_t   maxBottomToTopPhiH;

	Float_t maxBottomToTopRatioV;
	Int_t   maxBottomToTopPhiV;

} pLBlast;


void makeRecoTree(int run){

	AnitaVersion::set(4);

	AnitaDataset data(run,false,WaveCalType::kDefault,AnitaDataset::ANITA_ROOT_DATA,AnitaDataset::kRandomizePolarity);

	UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig(); 
	config->nmaxima = 3;
	//config->deconvolution_method = apd;
	config->combine_nantennas = 15;
	//config->correlator_theta_lowest = 42.0;

	UCorrelator::Analyzer *analyzer = new UCorrelator::Analyzer(config,true);
	UCorrelator::setAdaptiveFilterSpectrumAverageNSecs(10);
	UCorrelator::SineSubtractFilter::setUseCache(true);

	FilterStrategy *strategy = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");

	TString outname;
	outname.Form("recoTree_%d.root",run);
	TFile ofile(outname,"RECREATE");

	TTree * tree = new TTree("tree","Reconstructed pointed waveform tree");
	AnitaEventSummary * sum = new AnitaEventSummary; 

	hPol varsH;
	vPol varsV;	

	eventInfo info;
	pLBlast plb;

	tree->Branch("eventNumber",&info.eventNumber,"eventNumber/I");
	tree->Branch("run",&info.run,"run/I");
//	tree->Branch("isWais",&info.isWais,"isWais/I");
	tree->Branch("isRF",&info.isRF,"isRF/I");
	tree->Branch("isVarner",&info.isVarner,"isVarner/I");

	tree->Branch("maxBottomToTopRatioH",&plb.maxBottomToTopRatioH,"maxBottomToTopRatioH/F");
	tree->Branch("maxBottomToTopPhiH",&plb.maxBottomToTopPhiH,"maxBottomToTopPhiH/I");

	tree->Branch("thetaH0",&varsH.theta0,"thetaH0/F");	
	tree->Branch("thetaH1",&varsH.theta1,"thetaH1/F");	
	tree->Branch("thetaH2",&varsH.theta2,"thetaH2/F");	
	tree->Branch("phiH0",&varsH.phi0,"phiH0/F");	
	tree->Branch("phiH1",&varsH.phi1,"phiH1/F");	
	tree->Branch("phiH2",&varsH.phi2,"phiH2/F");	
	tree->Branch("mapPeakH0",&varsH.mapPeak0,"mapPeakH0/F");
	tree->Branch("mapPeakH1",&varsH.mapPeak1,"mapPeakH1/F");
	tree->Branch("mapPeakH2",&varsH.mapPeak2,"mapPeakH2/F");

	tree->Branch("coherent_hpeakH0",&varsH.coherent_hpeak0,"coherent_hpeakH0/F");
	tree->Branch("coherent_hpeakH1",&varsH.coherent_hpeak1,"coherent_hpeakH1/F");
	tree->Branch("coherent_hpeakH2",&varsH.coherent_hpeak2,"coherent_hpeakH2/F");
	tree->Branch("coherent_snrH0",&varsH.coherent_snr0,"coherent_snrH0/F");
	tree->Branch("coherent_snrH1",&varsH.coherent_snr1,"coherent_snrH1/F");
	tree->Branch("coherent_snrH2",&varsH.coherent_snr2,"coherent_snrH2/F");
	tree->Branch("coherent_filtered_hpeakH0",&varsH.coherent_filtered_hpeak0,"coherent_filtered_hpeakH0/F");
	tree->Branch("coherent_filtered_hpeakH1",&varsH.coherent_filtered_hpeak1,"coherent_filtered_hpeakH1/F");
	tree->Branch("coherent_filtered_hpeakH2",&varsH.coherent_filtered_hpeak2,"coherent_filtered_hpeakH2/F");
	tree->Branch("coherent_filtered_snrH0",&varsH.coherent_filtered_snr0,"coherent_filtered_snrH0/F");
	tree->Branch("coherent_filtered_snrH1",&varsH.coherent_filtered_snr1,"coherent_filtered_snrH1/F");
	tree->Branch("coherent_filtered_snrH2",&varsH.coherent_filtered_snr2,"coherent_filtered_snrH2/F");

	tree->Branch("deconv_hpeakH0",&varsH.deconv_hpeak0,"deconv_hpeakH0/F");
	tree->Branch("deconv_hpeakH1",&varsH.deconv_hpeak1,"deconv_hpeakH1/F");
	tree->Branch("deconv_hpeakH2",&varsH.deconv_hpeak2,"deconv_hpeakH2/F");
	tree->Branch("deconv_snrH0",&varsH.deconv_snr0,"deconv_snrH0/F");
	tree->Branch("deconv_snrH1",&varsH.deconv_snr1,"deconv_snrH1/F");
	tree->Branch("deconv_snrH2",&varsH.deconv_snr2,"deconv_snrH2/F");
	tree->Branch("deconv_filtered_hpeakH0",&varsH.deconv_filtered_hpeak0,"deconv_filtered_hpeakH0/F");
	tree->Branch("deconv_filtered_hpeakH1",&varsH.deconv_filtered_hpeak1,"deconv_filtered_hpeakH1/F");
	tree->Branch("deconv_filtered_hpeakH2",&varsH.deconv_filtered_hpeak2,"deconv_filtered_hpeakH2/F");
	tree->Branch("deconv_filtered_snrH0",&varsH.deconv_filtered_snr0,"deconv_filtered_snrH0/F");
	tree->Branch("deconv_filtered_snrH1",&varsH.deconv_filtered_snr1,"deconv_filtered_snrH1/F");
	tree->Branch("deconv_filtered_snrH2",&varsH.deconv_filtered_snr2,"deconv_filtered_snrH2/F");

	tree->Branch("maxBottomToTopRatioV",&plb.maxBottomToTopRatioV,"maxBottomToTopRatioV/F");
	tree->Branch("maxBottomToTopPhiV",&plb.maxBottomToTopPhiV,"maxBottomToTopPhiV/I");

	tree->Branch("thetaV0",&varsV.theta0,"thetaV0/F");	
	tree->Branch("thetaV1",&varsV.theta1,"thetaV1/F");	
	tree->Branch("thetaV2",&varsV.theta2,"thetaV2/F");	
	tree->Branch("phiV0",&varsV.phi0,"phiV0/F");	
	tree->Branch("phiV1",&varsV.phi1,"phiV1/F");	
	tree->Branch("phiV2",&varsV.phi2,"phiV2/F");	
	tree->Branch("mapPeakV0",&varsV.mapPeak0,"mapPeakV0/F");
	tree->Branch("mapPeakV1",&varsV.mapPeak1,"mapPeakV1/F");
	tree->Branch("mapPeakV2",&varsV.mapPeak2,"mapPeakV2/F");

	tree->Branch("coherent_hpeakV0",&varsV.coherent_hpeak0,"coherent_hpeakV0/F");
	tree->Branch("coherent_hpeakV1",&varsV.coherent_hpeak1,"coherent_hpeakV1/F");
	tree->Branch("coherent_hpeakV2",&varsV.coherent_hpeak2,"coherent_hpeakV2/F");
	tree->Branch("coherent_snrV0",&varsV.coherent_snr0,"coherent_snrV0/F");
	tree->Branch("coherent_snrV1",&varsV.coherent_snr1,"coherent_snrV1/F");
	tree->Branch("coherent_snrV2",&varsV.coherent_snr2,"coherent_snrV2/F");
	tree->Branch("coherent_filtered_hpeakV0",&varsV.coherent_filtered_hpeak0,"coherent_filtered_hpeakV0/F");
	tree->Branch("coherent_filtered_hpeakV1",&varsV.coherent_filtered_hpeak1,"coherent_filtered_hpeakV1/F");
	tree->Branch("coherent_filtered_hpeakV2",&varsV.coherent_filtered_hpeak2,"coherent_filtered_hpeakV2/F");
	tree->Branch("coherent_filtered_snrV0",&varsV.coherent_filtered_snr0,"coherent_filtered_snrV0/F");
	tree->Branch("coherent_filtered_snrV1",&varsV.coherent_filtered_snr1,"coherent_filtered_snrV1/F");
	tree->Branch("coherent_filtered_snrV2",&varsV.coherent_filtered_snr2,"coherent_filtered_snrV2/F");

	tree->Branch("deconv_hpeakV0",&varsV.deconv_hpeak0,"deconv_hpeakV0/F");
	tree->Branch("deconv_hpeakV1",&varsV.deconv_hpeak1,"deconv_hpeakV1/F");
	tree->Branch("deconv_hpeakV2",&varsV.deconv_hpeak2,"deconv_hpeakV2/F");
	tree->Branch("deconv_snrV0",&varsV.deconv_snr0,"deconv_snrV0/F");
	tree->Branch("deconv_snrV1",&varsV.deconv_snr1,"deconv_snrV1/F");
	tree->Branch("deconv_snrV2",&varsV.deconv_snr2,"deconv_snrV2/F");
	tree->Branch("deconv_filtered_hpeakV0",&varsV.deconv_filtered_hpeak0,"deconv_filtered_hpeakV0/F");
	tree->Branch("deconv_filtered_hpeakV1",&varsV.deconv_filtered_hpeak1,"deconv_filtered_hpeakV1/F");
	tree->Branch("deconv_filtered_hpeakV2",&varsV.deconv_filtered_hpeak2,"deconv_filtered_hpeakV2/F");
	tree->Branch("deconv_filtered_snrV0",&varsV.deconv_filtered_snr0,"deconv_filtered_snrV0/F");
	tree->Branch("deconv_filtered_snrV1",&varsV.deconv_filtered_snr1,"deconv_filtered_snrV1/F");
	tree->Branch("deconv_filtered_snrV2",&varsV.deconv_filtered_snr2,"deconv_filtered_snrV2/F");

	//variabl initialization//
	info.eventNumber               = -99999999;
	info.run                       = -99999999;
//	info.isWais                    = -99999999;
	info.isRF                      = -99999999;
	info.isVarner                  = -99999999;

	plb.maxBottomToTopRatioH       = -99999;
	plb.maxBottomToTopPhiH         = -99999;

	varsH.theta0                   = -99999;
	varsH.theta1                   = -99999;
	varsH.theta2                   = -99999;
	varsH.phi0                     = -99999;
	varsH.phi1                     = -99999;
	varsH.phi2                     = -99999;
	varsH.mapPeak0                 = -99999;
	varsH.mapPeak1                 = -99999;
	varsH.mapPeak2                 = -99999;

	varsH.coherent_hpeak0          = -99999;
	varsH.coherent_hpeak1          = -99999;
	varsH.coherent_hpeak2          = -99999;
	varsH.coherent_snr0            = -99999;
	varsH.coherent_snr1            = -99999;
	varsH.coherent_snr2            = -99999;
	varsH.coherent_filtered_hpeak0 = -99999;
	varsH.coherent_filtered_hpeak1 = -99999;
	varsH.coherent_filtered_hpeak2 = -99999;
	varsH.coherent_filtered_snr0   = -99999;
	varsH.coherent_filtered_snr1   = -99999;
	varsH.coherent_filtered_snr2   = -99999;

	varsH.deconv_hpeak0            = -99999;
	varsH.deconv_hpeak1            = -99999;
	varsH.deconv_hpeak2            = -99999;
	varsH.deconv_snr0              = -99999;
	varsH.deconv_snr1              = -99999;
	varsH.deconv_snr2              = -99999;
	varsH.deconv_filtered_hpeak0   = -99999;
	varsH.deconv_filtered_hpeak1   = -99999;
	varsH.deconv_filtered_hpeak2   = -99999;
	varsH.deconv_filtered_snr0     = -99999;
	varsH.deconv_filtered_snr1     = -99999;
	varsH.deconv_filtered_snr2     = -99999;

	plb.maxBottomToTopRatioV       = -99999;
	plb.maxBottomToTopPhiV         = -99999;

	varsV.theta0                   = -99999;
	varsV.theta1                   = -99999;
	varsV.theta2                   = -99999;
	varsV.phi0                     = -99999;
	varsV.phi1                     = -99999;
	varsV.phi2                     = -99999;
	varsV.mapPeak0                 = -99999;
	varsV.mapPeak1                 = -99999;
	varsV.mapPeak2                 = -99999;

	varsV.coherent_hpeak0          = -99999;
	varsV.coherent_hpeak1          = -99999;
	varsV.coherent_hpeak2          = -99999;
	varsV.coherent_snr0            = -99999;
	varsV.coherent_snr1            = -99999;
	varsV.coherent_snr2            = -99999;
	varsV.coherent_filtered_hpeak0 = -99999;
	varsV.coherent_filtered_hpeak1 = -99999;
	varsV.coherent_filtered_hpeak2 = -99999;
	varsV.coherent_filtered_snr0   = -99999;
	varsV.coherent_filtered_snr1   = -99999;
	varsV.coherent_filtered_snr2   = -99999;

	varsV.deconv_hpeak0            = -99999;
	varsV.deconv_hpeak1            = -99999;
	varsV.deconv_hpeak2            = -99999;
	varsV.deconv_snr0              = -99999;
	varsV.deconv_snr1              = -99999;
	varsV.deconv_snr2              = -99999;
	varsV.deconv_filtered_hpeak0   = -99999;
	varsV.deconv_filtered_hpeak1   = -99999;
	varsV.deconv_filtered_hpeak2   = -99999;
	varsV.deconv_filtered_snr0     = -99999;
	varsV.deconv_filtered_snr1     = -99999;
	varsV.deconv_filtered_snr2     = -99999;

	for (int i = 0; i < 100; i++){

		data.getEntry(i); 

		FilteredAnitaEvent ev(data.useful(), strategy, data.gps(), data.header()); 
		analyzer->analyze(&ev, sum);

		//hdr = data.header();
		//pat = data.gps(); 

		//if(-1*analyzer->getSummary()->peak[0][0].theta < 0.0) continue;//using physics convention for elevation angle
		//if(analyzer->getSummary()->flags.maxBottomToTopRatio[0] > 2.7) continue; 		
		//if(analyzer->getSummary()->flags.nSectorsWhereBottomExceedsTop > 26) continue;
		//if(analyzer->getSummary()->flags.meanPower[1+AnitaRing::kMiddleRing] / analyzer->getSummary()->flags.meanPower[1+AnitaRing::kMiddleRing] > 1.5) continue;
		//if(analyzer->getSummary()->flags.meanPower[1+AnitaRing::kBottomRing] / analyzer->getSummary()->flags.meanPower[1+AnitaRing::kMiddleRing] > 1.8) continue;
		//if(analyzer->getSummary()->flags.isPayloadBlast) continue;

		info.eventNumber               = analyzer->getSummary()->eventNumber;
		info.run                       = analyzer->getSummary()->run;

		//info.isWais                    = analyzer->getSummary()->flags.isWais;
		info.isRF                      = analyzer->getSummary()->flags.isRF;
		info.isVarner                  = analyzer->getSummary()->flags.isVarner;

		plb.maxBottomToTopRatioH       = analyzer->getSummary()->flags.maxBottomToTopRatio[0];
		plb.maxBottomToTopPhiH         = analyzer->getSummary()->flags.maxBottomToTopRatioSector[0];

		varsH.theta0                   = analyzer->getSummary()->peak[0][0].theta;
		varsH.theta1                   = analyzer->getSummary()->peak[0][1].theta;
		varsH.theta2                   = analyzer->getSummary()->peak[0][2].theta;
		varsH.phi0                     = analyzer->getSummary()->peak[0][0].phi;
		varsH.phi1                     = analyzer->getSummary()->peak[0][1].phi;
		varsH.phi2                     = analyzer->getSummary()->peak[0][2].phi;
		varsH.mapPeak0                 = analyzer->getSummary()->peak[0][0].value;
		varsH.mapPeak1                 = analyzer->getSummary()->peak[0][1].value;
		varsH.mapPeak2                 = analyzer->getSummary()->peak[0][2].value;

		varsH.coherent_hpeak0          = analyzer->getSummary()->coherent[0][0].peakHilbert;
		varsH.coherent_hpeak1          = analyzer->getSummary()->coherent[0][1].peakHilbert;
		varsH.coherent_hpeak2          = analyzer->getSummary()->coherent[0][2].peakHilbert;
		varsH.coherent_snr0            = analyzer->getSummary()->coherent[0][0].snr;
		varsH.coherent_snr1            = analyzer->getSummary()->coherent[0][1].snr;
		varsH.coherent_snr2            = analyzer->getSummary()->coherent[0][2].snr;
		varsH.coherent_filtered_hpeak0 = analyzer->getSummary()->coherent_filtered[0][0].peakHilbert;
		varsH.coherent_filtered_hpeak1 = analyzer->getSummary()->coherent_filtered[0][1].peakHilbert;
		varsH.coherent_filtered_hpeak2 = analyzer->getSummary()->coherent_filtered[0][2].peakHilbert;
		varsH.coherent_filtered_snr0   = analyzer->getSummary()->coherent_filtered[0][0].snr;
		varsH.coherent_filtered_snr1   = analyzer->getSummary()->coherent_filtered[0][1].snr;
		varsH.coherent_filtered_snr2   = analyzer->getSummary()->coherent_filtered[0][2].snr;
		
		varsH.deconv_hpeak0          = analyzer->getSummary()->deconvolved[0][0].peakHilbert;
		varsH.deconv_hpeak1          = analyzer->getSummary()->deconvolved[0][1].peakHilbert;
		varsH.deconv_hpeak2          = analyzer->getSummary()->deconvolved[0][2].peakHilbert;
		varsH.deconv_snr0            = analyzer->getSummary()->deconvolved[0][0].snr;
		varsH.deconv_snr1            = analyzer->getSummary()->deconvolved[0][1].snr;
		varsH.deconv_snr2            = analyzer->getSummary()->deconvolved[0][2].snr;
		varsH.deconv_filtered_hpeak0 = analyzer->getSummary()->deconvolved_filtered[0][0].peakHilbert;
		varsH.deconv_filtered_hpeak1 = analyzer->getSummary()->deconvolved_filtered[0][1].peakHilbert;
		varsH.deconv_filtered_hpeak2 = analyzer->getSummary()->deconvolved_filtered[0][2].peakHilbert;
		varsH.deconv_filtered_snr0   = analyzer->getSummary()->deconvolved_filtered[0][0].snr;
		varsH.deconv_filtered_snr1   = analyzer->getSummary()->deconvolved_filtered[0][1].snr;
		varsH.deconv_filtered_snr2   = analyzer->getSummary()->deconvolved_filtered[0][2].snr;

		plb.maxBottomToTopRatioV       = analyzer->getSummary()->flags.maxBottomToTopRatio[1];
		plb.maxBottomToTopPhiV         = analyzer->getSummary()->flags.maxBottomToTopRatioSector[1];

		varsV.theta0                   = analyzer->getSummary()->peak[1][0].theta;
		varsV.theta1                   = analyzer->getSummary()->peak[1][1].theta;
		varsV.theta2                   = analyzer->getSummary()->peak[1][2].theta;
		varsV.phi0                     = analyzer->getSummary()->peak[1][0].phi;
		varsV.phi1                     = analyzer->getSummary()->peak[1][1].phi;
		varsV.phi2                     = analyzer->getSummary()->peak[1][2].phi;
		varsV.mapPeak0                 = analyzer->getSummary()->peak[1][0].value;
		varsV.mapPeak1                 = analyzer->getSummary()->peak[1][1].value;
		varsV.mapPeak2                 = analyzer->getSummary()->peak[1][2].value;

		varsV.coherent_hpeak0          = analyzer->getSummary()->coherent[1][0].peakHilbert;
		varsV.coherent_hpeak1          = analyzer->getSummary()->coherent[1][1].peakHilbert;
		varsV.coherent_hpeak2          = analyzer->getSummary()->coherent[1][2].peakHilbert;
		varsV.coherent_snr0            = analyzer->getSummary()->coherent[1][0].snr;
		varsV.coherent_snr1            = analyzer->getSummary()->coherent[1][1].snr;
		varsV.coherent_snr2            = analyzer->getSummary()->coherent[1][2].snr;
		varsV.coherent_filtered_hpeak0 = analyzer->getSummary()->coherent_filtered[1][0].peakHilbert;
		varsV.coherent_filtered_hpeak1 = analyzer->getSummary()->coherent_filtered[1][1].peakHilbert;
		varsV.coherent_filtered_hpeak2 = analyzer->getSummary()->coherent_filtered[1][2].peakHilbert;
		varsV.coherent_filtered_snr0   = analyzer->getSummary()->coherent_filtered[1][0].snr;
		varsV.coherent_filtered_snr1   = analyzer->getSummary()->coherent_filtered[1][1].snr;
		varsV.coherent_filtered_snr2   = analyzer->getSummary()->coherent_filtered[1][2].snr;

		varsV.deconv_hpeak0          = analyzer->getSummary()->deconvolved[1][0].peakHilbert;
		varsV.deconv_hpeak1          = analyzer->getSummary()->deconvolved[1][1].peakHilbert;
		varsV.deconv_hpeak2          = analyzer->getSummary()->deconvolved[1][2].peakHilbert;
		varsV.deconv_snr0            = analyzer->getSummary()->deconvolved[1][0].snr;
		varsV.deconv_snr1            = analyzer->getSummary()->deconvolved[1][1].snr;
		varsV.deconv_snr2            = analyzer->getSummary()->deconvolved[1][2].snr;
		varsV.deconv_filtered_hpeak0 = analyzer->getSummary()->deconvolved_filtered[1][0].peakHilbert;
		varsV.deconv_filtered_hpeak1 = analyzer->getSummary()->deconvolved_filtered[1][1].peakHilbert;
		varsV.deconv_filtered_hpeak2 = analyzer->getSummary()->deconvolved_filtered[1][2].peakHilbert;
		varsV.deconv_filtered_snr0   = analyzer->getSummary()->deconvolved_filtered[1][0].snr;
		varsV.deconv_filtered_snr1   = analyzer->getSummary()->deconvolved_filtered[1][1].snr;
		varsV.deconv_filtered_snr2   = analyzer->getSummary()->deconvolved_filtered[1][2].snr;

		ofile.cd();
		tree->Fill();

	}

	ofile.cd();
	ofile.Write();

}


int main (int argc, char **argv){

	int run = atoi(argv[1]);
	
	makeRecoTree(run);	

}
