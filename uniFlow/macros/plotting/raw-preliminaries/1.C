for(Int_t k=0;k<15;k++)
{
    //cout << fv2_Phi_gap08_cent0->GetBinCenter(k+1)<<", ";
    v22Gap08_Phi_cc020[k] = fv2_Phi_gap08_cent0->GetBinContent(k+1);
    v22Gap08Err_Phi_cc020[k] = fv2_Phi_gap08_cent0->GetBinError(k+1);
    v22Gap08Sys_Phi_cc020[k] = sys_Phi * fv2_Phi_gap08_cent0->GetBinContent(k+1);
    
}

TGraphErrors *figv2_Phi_gap08_cent0 = new TGraphErrors(18,pt_Phi,v22Gap08_Phi_cc020,ptBin, v22Gap08Err_Phi_cc020);
figv2_Phi_gap08_cent0->SetMarkerStyle(34);
figv2_Phi_gap08_cent0->SetMarkerColor(kGreen+3);
figv2_Phi_gap08_cent0->SetMarkerSize(1.3);
figv2_Phi_gap08_cent0->SetLineColor(kGreen+3);
figv2_Phi_gap08_cent0->SetLineWidth(2);
figv2_Phi_gap08_cent0->SetLineStyle(1);

TGraphErrors *gv2_Phi_gap08_cent0 = new TGraphErrors(18,pt_Phi,v22Gap08_Phi_cc020,ptErr, v22Gap08Err_Phi_cc020);
gv2_Phi_gap08_cent0->SetMarkerStyle(34);
gv2_Phi_gap08_cent0->SetMarkerColor(kGreen+3);
gv2_Phi_gap08_cent0->SetMarkerSize(1.3);
gv2_Phi_gap08_cent0->SetLineColor(kGreen+3);
gv2_Phi_gap08_cent0->SetLineWidth(2);
gv2_Phi_gap08_cent0->SetLineStyle(1);
gv2_Phi_gap08_cent0->SetFillStyle(1001);
gv2_Phi_gap08_cent0->SetFillColor(kGreen-8);


gv2_Phi_gap08_cent0->Draw("2psame");
figv2_Phi_gap08_cent0->Draw("pZsame");
