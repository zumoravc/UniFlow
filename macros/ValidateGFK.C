void ValidateGFK()
{
	//TFile* fInput = new TFile("~/NBI/Flow/results/V0s/GFK-diff/merge/AnalysisResults.root","READ");
	TFile* fInput = new TFile("~/NBI/Flow/flow/AnalysisResults.root","READ");
	fInput->cd("flowPID_FB768");
	//gDirectory->ls();

	TList* lPID = (TList*) gDirectory->Get("PID_flowPID_FB768");
	TList* lTracks = (TList*) gDirectory->Get("Tracks_flowPID_FB768");
	TList* lV0s = (TList*) gDirectory->Get("V0s_flowPID_FB768");
	//lPID->ls();

	//TProfile* profOld = (TProfile*) lTracks->FindObject("fTracksRefTwo_n3_Gap10_sample0");
	//TProfile* profOld = (TProfile*) lPID->FindObject("fPionsDiffTwoPos_n2_Gap-10_Cent1_sample0");
	//TProfile* profOld2 = (TProfile*) lPID->FindObject("fKaonsDiffTwoNeg_n3_Gap00_Cent3_sample0");
	//TProfile* profOld3 = (TProfile*) lPID->FindObject("fProtonsDiffTwoPos_n4_Gap10_Cent4_sample0");
	//TProfile* profOld4 = (TProfile*) lPID->FindObject("fProtonsDiffTwoPos_n2_Gap00_Cent2_sample3");
	//TProfile* profOld5 = (TProfile*) lPID->FindObject("fProtonsDiffTwoPos_n2_Gap00_Cent2_sample4");
	
	TProfile* profOld = (TProfile*) lTracks->FindObject("fTracksDiffTwoPos_n2_Gap-10_Cent1_sample0");
	TProfile* profOld2 = (TProfile*) lTracks->FindObject("fTracksDiffTwoPos_n2_Gap10_Cent3_sample1");
	TProfile* profOld3 = (TProfile*) lTracks->FindObject("fTracksDiffTwoPos_n3_Gap-10_Cent3_sample2");
	TProfile* profOld4 = (TProfile*) lV0s->FindObject("fV0sDiffTwoPos_Lambda_n2_Gap-10_Cent3_sample3");
	TProfile* profOld5 = (TProfile*) lV0s->FindObject("fV0sDiffTwoPos_Lambda_n2_Gap-10_Cent3_sample4");
	
	//TProfile* profGFK = (TProfile*) lPID->FindObject("fPion_n22_Pos_gap-10_cent1_number0");
	//TProfile* profGFK2 = (TProfile*) lPID->FindObject("fKaon_n32_Neg_gap00_cent3_number0");
	//TProfile* profGFK3 = (TProfile*) lPID->FindObject("fProton_n42_Pos_gap10_cent4_number0");
	
	TProfile* profGFK = (TProfile*) lTracks->FindObject("fTracks_n22_Pos_gap-10_cent1_number0");
	TProfile* profGFK2 = (TProfile*) lTracks->FindObject("fTracks_n22_Pos_gap10_cent3_number1");
	TProfile* profGFK3 = (TProfile*) lTracks->FindObject("fTracks_n32_Pos_gap-10_cent3_number2");
	//TProfile* profGFK4 = (TProfile*) lPID->FindObject("fd22Proton_gap00_cent1_number3");
	//TProfile* profGFK5 = (TProfile*) lPID->FindObject("fd22Proton_gap00_cent1_number4");

	TCanvas* cCan = new TCanvas();
	cCan->Divide(2,3);
	cCan->cd(1);
	profOld->Draw("colz");
	cCan->cd(2);
	profGFK->Draw("colz");
	cCan->cd(3);
	profOld2->Draw("colz");
	cCan->cd(4);
	profGFK2->Draw("colz");
	//profOld3->Draw("colz");
	cCan->cd(5);
	profOld3->Draw("colz");
	cCan->cd(6);
	//profOld5->Draw("colz");
	profGFK3->Draw("colz");
	
}