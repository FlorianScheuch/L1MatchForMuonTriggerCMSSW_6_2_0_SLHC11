{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Mon Sep  8 16:32:51 2014) by ROOT version5.34/07
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",258,102,1483,944);
   Canvas_1->Range(-0.125,-0.8693953,6.125,4.814257);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetLogy();
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   
   TH1D *NumberOf L1 muon matches = new TH1D("NumberOf L1 muon matches","Number of L1 muon matches",20,-0.5,19.5);
   NumberOf L1 muon matches->SetBinContent(2,9297);
   NumberOf L1 muon matches->SetBinContent(3,281);
   NumberOf L1 muon matches->SetBinContent(4,9);
   NumberOf L1 muon matches->SetBinContent(5,1);
   NumberOf L1 muon matches->SetEntries(9588);
   NumberOf L1 muon matches->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   NumberOf L1 muon matches->SetLineColor(ci);
   NumberOf L1 muon matches->SetLineWidth(4);
   NumberOf L1 muon matches->GetXaxis()->SetTitle("# of matched L1 muons to one GEN muon");
   NumberOf L1 muon matches->GetXaxis()->SetRange(2,6);
   NumberOf L1 muon matches->GetXaxis()->SetNdivisions(5);
   NumberOf L1 muon matches->GetXaxis()->SetLabelFont(42);
   NumberOf L1 muon matches->GetXaxis()->SetLabelSize(0.05);
   NumberOf L1 muon matches->GetXaxis()->SetTitleSize(0.05);
   NumberOf L1 muon matches->GetXaxis()->SetTitleFont(42);
   NumberOf L1 muon matches->GetYaxis()->SetTitle("#");
   NumberOf L1 muon matches->GetYaxis()->SetLabelFont(42);
   NumberOf L1 muon matches->GetYaxis()->SetLabelSize(0.06);
   NumberOf L1 muon matches->GetYaxis()->SetTitleSize(0.06);
   NumberOf L1 muon matches->GetYaxis()->SetTitleOffset(0.83);
   NumberOf L1 muon matches->GetYaxis()->SetTitleFont(42);
   NumberOf L1 muon matches->GetZaxis()->SetLabelFont(42);
   NumberOf L1 muon matches->GetZaxis()->SetLabelSize(0.035);
   NumberOf L1 muon matches->GetZaxis()->SetTitleSize(0.035);
   NumberOf L1 muon matches->GetZaxis()->SetTitleFont(42);
   NumberOf L1 muon matches->Draw("");
   
   TPaveText *pt = new TPaveText(0.2919952,0.94,0.7582426,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *text = pt->AddText("Number of L1 muon matches");
   pt->Draw();
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
