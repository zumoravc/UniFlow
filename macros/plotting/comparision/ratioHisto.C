void ratioHisto()
{
  TFile* fileIn = TFile::Open("/home/alidock/ana/output/test/Processed.root","READ");

  Double_t xCross[10] = {0,5,10,20,30,40,50,60,70,80};
  // RUN2 data
  Double_t v22Gap10Run2[10]={0.0283859, 0.0456604, 0.0655068, 0.0870721, 0.099105, 0.104143, 0.10286, 0.0974591, 0.0888104};
  Double_t v22Gap10Run2CombErr[10]={0.000711585, 0.000940029, 0.00105135, 0.00137951, 0.00158637, 0.00172523, 0.00187963, 0.00236759, 0.0045766};

  TH1D *data_v22Gap10 = (TH1D*) fileIn->Get("Refs_hFlow2_harm2_gap10");
  if(!data_v22Gap10) {printf("OOPS \n");}; return;}
  //TH1D *published_v22 = new TH1D("published_v22", "published", 9, xCross);
  TH1D *published_v22 = (TH1D*)data_v22Gap10->Clone("published_v22");
  for(int i = 1; i < 10; i++)
  {
    published_v22->SetBinContent(i,v22Gap10Run2[i-1]);
    published_v22->SetBinError(i,v22Gap10Run2CombErr[i-1]);
  }
  TH1D *ratio = (TH1D*)data_v22Gap10->Clone("ratio");
  ratio->Divide(published_v22);
  ratio->Draw();


}
