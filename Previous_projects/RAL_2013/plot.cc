//ROOT HEADERS
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TLorentzVector.h>


//MY OBJECT HEADERS
#include "TwoBody.hh"

//STL HEADERS
#include <iostream>
#include <math.h>
#include <vector>
#include <string>

using std::cout;
using std::cerr;
using std::vector;
using std::string;

using eol::TwoBody;
using eol::GenParticle;

void mypause()
{ 
  std::cout<<"Press [Enter] to continue . . .";
  std::cin.get();
} 


int main(int argc, char **argv)
{

  if (argc != 2) {
    cerr << "Wrong number of arguments\n";
    cerr << "eg. Should be\n./plot <myinput.root>\n";
    cout << argc << "\n";
    cout << argv[1] << "\n";
    return -1;
  }

  int arraySize;

  {
    //create a large character string on the stack
    //it is in a local space so will dissapear
    //once we know the size fo the char string
    char buffer[1000];
    arraySize = sprintf(buffer,argv[1]);
  }

  char fileName[arraySize];
  sprintf(fileName,argv[1]);

  //open ROOT file
  TFile* f = new TFile(fileName);

  //BEWARE, TApplication changes argc and i think argv
  //make sure you have finished with argc and argv before using this
  TApplication graphicsPlease("graphicsPlease", &argc, argv); 
  
  //Assign a TTree pointer to tree
  TTree *evtTree = (TTree*)f->Get("evtTree");
  
  //Lets see what is inside the tree
  evtTree->Print();
 
  //first create a canvas
  TCanvas* massCan = new TCanvas("massCan","tb mass distro");
  massCan->cd(); //move into canvas
  evtTree->Draw("m"); //draw Two Body mass 


  //Define Branch, TwoBody and  GenParticle.
  TwoBody* tbEvent = new TwoBody();
  GenParticle* xEvent = new GenParticle();
  GenParticle* missEvent = new GenParticle();
  GenParticle* inpEvent = new GenParticle(); 
  TBranch* branch = evtTree->GetBranch("tbEvent");
  TBranch* xbranch = evtTree->GetBranch("xEvent");
  TBranch* missbranch = evtTree->GetBranch("missEvent");
  TBranch* inpbranch = evtTree->GetBranch("inpEvent");
  const GenParticle* gen1 ;
  const GenParticle* gen2 ;  
  branch->SetAddress(&tbEvent); //set loaded value to tbEvent pointer. Confusing as I have used the same name
  xbranch->SetAddress(&xEvent);
  missbranch->SetAddress(&missEvent);
  inpbranch->SetAddress(&inpEvent);
  
  
  //Set Up Lorentz Vectors
  TLorentzVector* gen1_lz = new TLorentzVector;
  TLorentzVector* gen2_lz = new TLorentzVector;
 TLorentzVector* xEvent_lz = new TLorentzVector;  

  //Set Parameters for Histos   
  unsigned int nBins = 500; // Number of Bins for general histos
  unsigned int nBins_sp = 100; // Number of Bins for specific event histos
  unsigned int mMax = 2; //MaxMass
  unsigned int mMin = 0; //MinMass
  unsigned int xmax = 1000; //Max Energy/TwoBody Mass
  unsigned int xmin = 0; //Min Energy/TwoBody Mass

  //number of events in tree
  const unsigned int numEvents = branch->GetEntries();
   
  bool fit;
  cout << "\nDo you to fit the histograms? (1 or 0)\n";
  std::cin >> fit; 

  bool specificEvent= 0;      // Used to decide wether to plot specific events. 
  cout << "\nDo you want to plot specific events (1 or 0) [E.g. electron-electron, electron-muon ...]\n";
  std::cin >> specificEvent;  

  //Declare Histos

  //TwoBody Histogams
  TH1F histo_tb_m("histo_tb_m", "Two Body Mass Histo", nBins, xmin, xmax);
  TH1F histo_tb_E("histo_tb_E", "Two Body Energy Histo", nBins, xmin, xmax);
  TH1F histo_tb_pt("histo_tb_pt","Two Body pt Histo", nBins, xmin, xmax);
  TH1F histo_tb_ptot("histo_tb_ptot","Two Body total p Histo", nBins, xmin, xmax);
  TH1F histo_tb_th("histo_tb_th", "Two Body theta Histo",40, 0, 3.5);
  TH1F histo_tb_costh("histo_tb_costh", "Two Body cos(theta) Histo", 40, -1 , 1);

  //Input Histograms
  TH1F histo_input_E("histo_input_E", "Input Energy Histo", nBins, 0, 4000);
  TH1F histo_input_ptot("histo_input_ptot","Input total p Histo", nBins, 0, 4000);
  
  //Lepton Histograms
  TH1F histo_leptot_m("histo_lt_m","Total Lepton Mass Histo",nBins,mMin-0.2,mMax); //Histogram for both lepton 1 and lepton 2
  TH1F histo_mSum("histo_mSum","Sum of Masses of lep1 and lep2", nBins,mMin-0.2,2*mMax); //Histogram of sum of masses
 
  //Missing Energy Histograms
  TH1F histo_miss_m("histo_miss_m", "Missing Energy Mass Histo", nBins, xmin, xmax);
  TH1F histo_miss_E("histo_miss_E", "Missing Energy Energy Histo", nBins, xmin, xmax);
  TH1F histo_miss_pt("histo_miss_pt","Missing Energy pt Histo", nBins, xmin, xmax);

  // X Frame events
  TH1F histo_xframe_pt1("histo_xframe_pt1", "pt of gen1 in x frame", 1000, 0, 300);
  TH1F histo_xframe_pt2("histo_xframe_pt2", "pt of gen2 in x frame", 1000, 0, 300);
  TH1F histo_xframe_theta("histo_xframe_theta", "Theta between particles in x frame", 1000, 0, 3.25);
  TH1F histo_xframe_costh("histo_xframe_costh", "Cos(theta) in x frame", 1000, -1.1, 1.1);

  // Specific Events -> They have three histos for twobody mass, energy and pt
  // Electron-Electron Events
  TH1F histo_ee_m("hisito_ee_m", "Electron-Electron Mass Histo", nBins_sp, xmin, xmax);
  TH1F histo_ee_E("hisito_ee_E", "Electron-Electron Energy Histo", nBins_sp, xmin, xmax);
  TH1F histo_ee_pt("hisito_ee_pt", "Electron-Electron pt Histo", nBins_sp, xmin, xmax);
  TH1F histo_ee_th("histo_ee_th", "Electron-Electron theta Histo",40, 0, 3.5);
  TH1F histo_ee_costh("histo_ee_costh", "Electron-Electron cos(theta) Histo",40,-1 ,1);

  // Electron-Mu Events
  TH1F histo_emu_m("hisito_emu_m", "Electron-Muon Mass Histo", nBins_sp, xmin, xmax);
  TH1F histo_emu_E("hisito_emu_E", "Electron-Muon Energy Histo", nBins_sp, xmin, xmax);
  TH1F histo_emu_pt("hisito_emu_pt", "Electron-Muon pt Histo", nBins_sp, xmin, xmax);
  TH1F histo_emu_th("histo_emu_th", "Electron-Muon theta Histo",40, 0, 3.5);
  TH1F histo_emu_costh("histo_emu_costh", "Electron-Muon cos(theta) Histo",40,-1 ,1);

  // Electron-Tau Events
  TH1F histo_etau_m("hisito_etau_m", "Electron-Tau Mass Histo", nBins_sp, xmin, xmax);
  TH1F histo_etau_E("hisito_etau_E", "Electron-Tau Energy Histo", nBins_sp, xmin, xmax);
  TH1F histo_etau_pt("hisito_etau_pt", "Electron-Tau pt Histo", nBins_sp, xmin, xmax); 
  TH1F histo_etau_th("histo_etau_th", "Electron-Tau theta Histo",40, 0, 3.5);
  TH1F histo_etau_costh("histo_etau_costh", "Electron-Tau cos(theta) Histo",40,-1 ,1);

  // Muon-Muon Events
  TH1F histo_mumu_m("hisito_mumu_m", "Muon-Muon Mass Histo", nBins_sp, xmin, xmax);
  TH1F histo_mumu_E("hisito_mumu_E", "Muon-Muon Energy Histo", nBins_sp, xmin, xmax);
  TH1F histo_mumu_pt("hisito_mumu_pt", "Muon-Muon pt Histo", nBins_sp, xmin, xmax); 
  TH1F histo_mumu_th("histo_mumu_th", "Muon-Muon theta Histo",40, 0, 3.5);
  TH1F histo_mumu_costh("histo_mumu_costh", "Muon-Muon cos(theta) Histo",40,-1 ,1);

  // Muon-Tau Events
  TH1F histo_mutau_m("hisito_mutau_m", "Muon-Tau Mass Histo", nBins_sp, xmin, xmax);
  TH1F histo_mutau_E("hisito_mutau_E", "Muon-Tau Energy Histo", nBins_sp, xmin, xmax);
  TH1F histo_mutau_pt("hisito_mutau_pt", "Muon-Tau pt Histo", nBins_sp, xmin, xmax); 
  TH1F histo_mutau_th("histo_mutau_th", "Muon-Tau theta Histo",40, 0, 3.5);
  TH1F histo_mutau_costh("histo_mutau_costh", "Muon-Tau cos(theta) Histo",40,-1 ,1);

  // TauTau Events
  TH1F histo_tautau_m("hisito_tautau_m", "Tau-Tau Mass Histo", nBins_sp, xmin, xmax);
  TH1F histo_tautau_E("hisito_tautau_E", "Tau-Tau Energy Histo", nBins_sp, xmin, xmax);
  TH1F histo_tautau_pt("hisito_tautau_pt", "Tau-Tau pt Histo", nBins_sp, xmin, xmax); 
  TH1F histo_tautau_th("histo_tautau_th", "Tau-Tau theta Histo",40, 0, 3.5);
  TH1F histo_tautau_costh("histo_tautau_costh", "Tau-Tau cos(theta) Histo",40,-1 ,1);

   
  //loop over tree and fill histo 
  for (unsigned int i = 0; i < numEvents; ++i){
    branch->GetEntry(i); // set tree object tbEvent for each event i
    xbranch->GetEntry(i);
    missbranch->GetEntry(i);
    inpbranch->GetEntry(i);
    gen1 = tbEvent->getGen1();  // Retrieve Lepton 1 
    gen2 = tbEvent->getGen2();  // Retrieve Lepton 2 
    //Find the sum of masses
    double mass_sum;
    mass_sum = (gen1->getm())+(gen2->getm());
    
    double px_1 = (gen1->getpx());
    double py_1 = (gen1->getpy());
    double pz_1 = (gen1->getpz());
    double px_2 = (gen2->getpx());
    double py_2 = (gen2->getpy());
    double pz_2 = (gen2->getpz());
    double px_x = (xEvent->getpx());
    double py_x = (xEvent->getpy());
    double pz_x = (xEvent->getpz());    
    double theta;
    double cos_th;
    
    cos_th = (((px_1*px_2)+(py_1*py_2)+(pz_1*pz_2))/((sqrt((px_1*px_1)+(py_1*py_1)+(pz_1*pz_1)))*(sqrt((px_2*px_2)+(py_2*py_2)+(pz_2*pz_2))))); 
    theta = acos(cos_th);
    // cout << "\n Theta = " << theta <<"    \t" << i;     

    double ptot = sqrt((tbEvent->getpx())*(tbEvent->getpx())+(tbEvent->getpy())*(tbEvent->getpy())+(tbEvent->getpz())*(tbEvent->getpz()));
    double ptot_inp = sqrt((inpEvent->getpx())*(inpEvent->getpx())+(inpEvent->getpy())*(inpEvent->getpy())+(inpEvent->getpz())*(inpEvent->getpz()));

    //Fill TwoBody Histos
    histo_tb_m.Fill(tbEvent->getm());    
    histo_tb_E.Fill(tbEvent->getE()); 
    histo_tb_pt.Fill(tbEvent->getpt());
    histo_tb_ptot.Fill(ptot);
    histo_tb_th.Fill(theta);
    histo_tb_costh.Fill(cos_th);

    //Fill input histos
    histo_input_E.Fill(abs(inpEvent->getE()));
    histo_input_ptot.Fill(ptot_inp);
    
    //Fill General Lepton Histos
    histo_leptot_m.Fill(gen1->getm());  
    histo_leptot_m.Fill(gen2->getm());
    histo_mSum.Fill(mass_sum);
        
    //Fill missing energy histos
    histo_miss_m.Fill(missEvent->getm());
    histo_miss_E.Fill(missEvent->getE());
    histo_miss_pt.Fill(missEvent->getpt());
        
    //Fill Lorentz vectors with (px,py,pz,E)
    gen1_lz->SetPxPyPzE(px_1,py_1,pz_1,gen1->getE());
    gen2_lz->SetPxPyPzE(px_2,py_2,pz_2,gen2->getE()); 
    xEvent_lz->SetPxPyPzE(px_x,py_x,pz_x,xEvent->getE()); 

    // Boost into X particle frame
    gen1_lz->Boost(-xEvent_lz->BoostVector()); 
    gen2_lz->Boost(-xEvent_lz->BoostVector());  
    xEvent_lz->Boost(-xEvent_lz->BoostVector());
     
    double px_lz_1 = (gen1_lz->Px());
    double py_lz_1 = (gen1_lz->Py());
    double pz_lz_1 = (gen1_lz->Pz());
    double px_lz_2 = (gen2_lz->Px());
    double py_lz_2 = (gen2_lz->Py());
    double pz_lz_2 = (gen2_lz->Pz());
    double px_lz_x = (xEvent_lz->Px());
    double py_lz_x = (xEvent_lz->Py());
    double pz_lz_x = (xEvent_lz->Pz());    
    double theta_lz;
    double cos_th_lz;

    cos_th_lz = (((px_lz_1*px_lz_2)+(py_lz_1*py_lz_2)+(pz_lz_1*pz_lz_2)) // runs over line
             / ((sqrt((px_lz_1*px_lz_1)+(py_lz_1*py_lz_1)+(pz_lz_1*pz_lz_1)))*(sqrt((px_lz_2*px_lz_2)+(py_lz_2*py_lz_2)+(pz_lz_2*pz_lz_2))))); 
    theta_lz = acos(cos_th_lz);     

    //Fill X frame
    histo_xframe_pt1.Fill(sqrt(px_lz_1*px_lz_1 + py_lz_1*py_lz_1));
    histo_xframe_pt2.Fill(sqrt(px_lz_2*px_lz_2 + py_lz_2*py_lz_2));
    histo_xframe_theta.Fill(theta_lz);
    histo_xframe_costh.Fill(cos_th_lz);
    // cout << "\t" << test_px;

   if (specificEvent == 1){
    //Fill Lepton Specific Events
    if (mass_sum == 2*0.000511) { //if two electrons fill ee histos
     histo_ee_m.Fill(tbEvent->getm());   
     histo_ee_E.Fill(tbEvent->getE());
     histo_ee_pt.Fill(tbEvent->getpt());
     histo_ee_th.Fill(theta);
     histo_ee_costh.Fill(cos_th);
    }   
    if (mass_sum == (0.000511+0.1056)) { //if electron-muon fill emu histos
     histo_emu_m.Fill(tbEvent->getm());   
     histo_emu_E.Fill(tbEvent->getE());
     histo_emu_pt.Fill(tbEvent->getpt());
     histo_emu_th.Fill(theta);
     histo_emu_costh.Fill(cos_th);
    }   
    if (mass_sum == (0.000511+1.777)) { //if electron-tau fill etau histos
     histo_etau_m.Fill(tbEvent->getm());   
     histo_etau_E.Fill(tbEvent->getE());
     histo_etau_pt.Fill(tbEvent->getpt());
     histo_etau_th.Fill(theta);
     histo_etau_costh.Fill(cos_th);
    }   
    if (mass_sum == (2*0.1056)) { //if two muons fill mumu histos
     histo_mumu_m.Fill(tbEvent->getm());   
     histo_mumu_E.Fill(tbEvent->getE());
     histo_mumu_pt.Fill(tbEvent->getpt());
     histo_mumu_th.Fill(theta);
     histo_mumu_costh.Fill(cos_th);
    }   
    if (mass_sum ==(0.1056+1.777)) { //if muon-tau fill mutau histos
     histo_mutau_m.Fill(tbEvent->getm());   
     histo_mutau_E.Fill(tbEvent->getE());
     histo_mutau_pt.Fill(tbEvent->getpt());
     histo_mutau_th.Fill(theta);
     histo_mutau_costh.Fill(cos_th);
    }   
    if (mass_sum == (2*1.777)) { //if two taus fill tautau histos
     histo_tautau_m.Fill(tbEvent->getm());   
     histo_tautau_E.Fill(tbEvent->getE());
     histo_tautau_pt.Fill(tbEvent->getpt());
     histo_tautau_th.Fill(theta);
     histo_tautau_costh.Fill(cos_th);
    }    
   }
  }
  
   
  if(fit == 1) {  

    //Define Vectors to store parameter data
    vector<double> a_m;
    vector<double> b_m;
    vector<double> n_m;
    vector<double> a_E;
    vector<double> n_E;
    vector<double> b_E;
    vector<double> a_pt;
    vector<double> n_pt;
    vector<double> b_pt;

   //Create Fitting Functions
   // f1=A*x*exp(-x/B)
   TF1 *f1 = new TF1("f1","[0]*(x)*exp(-x/[1])",0,1000);
     f1->SetParameters(999,999);  
   // f2=A*(x^2)*exp(-x/B)
   TF1 *f2 = new TF1("f2","[0]*(x^2)*exp(-x/[1])",0,1000);
     f2->SetParameters(999,999);
   // f3=A*(x^3)*exp(-x/B)
   TF1 *f3 = new TF1("f3","[0]*(x^3)*exp(-x/[1])",0,1000);
     f3->SetParameters(999,999);
   // f4=A*(x^4)*exp(-x/B)
   TF1 *f4 = new TF1("f4","[0]*(x^4)*exp(-x/[1])",0,1000);
     f4->SetParameters(999,999);
   // fn=A*(x^n)*exp(-x/B)
   TF1 *fn = new TF1("fn","[0]*(x^[2])*exp(-x/[1])",0,1000);
     fn->SetParameters(999,999,999);
   TF1 *fsin = new TF1("fsin","[0]*sin(x)",0, 3.5);
     fsin->SetParameter(0,999);
    
  
   TF1 *f_m = new TF1("f_m","f2",0,1000); 
     f_m->SetParNames("Const","Expo Divide","n, (x^n)");
     f_m->SetParameters(999,999);
   TF1 *f_E = new TF1("f_E","f3",0,1000); 
     f_E->SetParNames("Const","Expo Divide","n, (x^n)");
     f_E->SetParameters(999,999);
   TF1 *f_pt = new TF1("f_pt","f1",0,1000); 
     f_pt->SetParNames("Const","Expo Divide","n, (x^n)");
     f_pt->SetParameters(999,999);
 

   // Apply Fittings
   // Two Body event fitting
    cout<< "\n \n fit to histo_tb_m \n";      // Fit to mass histo
      histo_tb_m.Fit("f_m");
       a_m.push_back(f_m->GetParameter(0));     // Storing parameter a for the two body mass.
       b_m.push_back(f_m->GetParameter(1));     // Scoring parameter b for the two body mass.
       n_m.push_back(f_m->GetParameter(2));     // Scoring parameter n for the two body mass.
    cout<< "\n\n fit to histo_tb_E \n";       // Fit to Energy histo
      histo_tb_E.Fit("f_E");    
       a_E.push_back(f_E->GetParameter(0));
       b_E.push_back(f_E->GetParameter(1));
       n_E.push_back(f_E->GetParameter(2));
    cout<< "\n\n fit to histo_tb_pt \n";      // Fit to pt histo
      histo_tb_pt.Fit("f_pt");    
       a_pt.push_back(f_pt->GetParameter(0));
       b_pt.push_back(f_pt->GetParameter(1)); 
       n_pt.push_back(f_pt->GetParameter(2));
       


       cout << "\n\nParameter Fitting Data";
       cout <<"\nj \ta_m \t\tb_m \t\tn_m \t\ta_E \t\tb_E \t\tn_E  \t\ta_pt \t\tb_pt \t\tn_pt";
      for( int j =0; j<3; ++j) {
        cout << "\n" << j << "\t" << a_m[j] << "    \t" << b_m[j] << "   \t" << n_m[j]<< "\t\t" ;
        cout << a_E[j] <<"\t" << b_E[j] << "\t\t" << n_E[j]  << "\t\t" ;
        cout << a_pt[j] << "   \t" << b_pt[j] <<  "\t\t" << n_pt[j];
      }
  }
  else {cout << "\nFitting Not Done";}


   massCan->cd(); //move into canvas
   evtTree->Draw("m"); //draw Two Body mass   


  //Canvas containing TwoBody Events
  TCanvas* tbCan = new TCanvas("tbCan","tb distro total");
  tbCan->Divide(1,3);
  // Plot histos
  tbCan->cd(1);
    histo_tb_m.Draw("");
  tbCan->cd(2);
    histo_tb_E.Draw("");
  tbCan->cd(3);
    histo_tb_pt.Draw("");

  //Canvas for TwoBody total momentum
  TCanvas* pTotCan = new TCanvas("pTotCan", "total momentum distro");
  pTotCan->cd();
     //Next Three Lines are to normalise total plot
      //Double_t integral = histo_tb_ptot.Integral();
      //histo_tb_ptot.Scale(1.0/integral); 
      //cout << "\n integral:   " << integral; 
    histo_tb_ptot.Draw("");

  TCanvas* tbThetaCan = new TCanvas("tbThetaCan","tb theta distro total");
  tbThetaCan->Divide(1,2);
  tbThetaCan->cd(1);
  histo_tb_th.Draw("");
  tbThetaCan->cd(2);
  histo_tb_costh.Draw("");

   //Input canvas
  TCanvas* inpCan = new TCanvas("inpCan", "Input particles distro");
  inpCan->Divide(1,2);
  inpCan->cd(1);
  histo_input_E.Draw("");
  inpCan->cd(2); 
  histo_input_ptot.Draw("");
  
  //Missing Energy canvas
  TCanvas* missCan = new TCanvas("missCan", "Missing Energy distro");
  missCan->Divide(1,3);
  missCan->cd(1);
  histo_miss_m.Draw("");
  missCan->cd(2);
  histo_miss_E.Draw("");
  missCan->cd(3);
  histo_miss_pt.Draw("");

  //X frame canvas
  TCanvas* xframeCan = new TCanvas("xFrameCan", "X frame distro");
  xframeCan->Divide(1,4);
  xframeCan->cd(1);
  histo_xframe_pt1.Draw("");
  xframeCan->cd(2);
  histo_xframe_pt2.Draw("");
  xframeCan->cd(3);
  histo_xframe_theta.Draw("");
  xframeCan->cd(4);
  histo_xframe_costh.Draw("");

  //Analyse GenParticle 1 and GenParticle 2
  //New Canvas Divide in (1,2) Plot Lepton 1 charge and Lepton 2 charge
  TCanvas* lepCan = new TCanvas("lepCan", "lepton distro");
  lepCan-> Divide(2,3);
  //Plot
  lepCan-> cd(1);
  evtTree->Draw("gen1.charge");
  lepCan-> cd(2);
  evtTree->Draw("gen2.charge");
  lepCan->cd(3);
  evtTree->Draw("gen1.m");
  lepCan->cd(4);
  evtTree->Draw("gen2.m");
  lepCan->cd(5);
  evtTree->Draw("gen1.pt");
  lepCan->cd(6);
  evtTree->Draw("gen2.pt");
  
  //New Canvas and Draw total lepton mass
  TCanvas* lepTotCan = new TCanvas("lepTotCan", "Total Lepton Disto");
  lepTotCan->Divide(1,2);
  lepTotCan->cd(1);
  histo_leptot_m.Draw("");
  lepTotCan->cd(2);
  histo_mSum.Draw("");

  
  if (specificEvent == 1) {
  //Canvas for ee
  TCanvas* eeCan = new TCanvas("eeCan", "Electron-Electron Events");
   eeCan->Divide(1,3);
   eeCan->cd(1);
     histo_ee_m.Draw("");
   eeCan->cd(2);
     histo_ee_E.Draw("");
   eeCan->cd(3);
     histo_ee_pt.Draw("");
  //Canvas for eMu
  TCanvas* emuCan = new TCanvas("emuCan", "Electron-Muon Events");
   emuCan->Divide(1,3);
   emuCan->cd(1);
     histo_emu_m.Draw("");
   emuCan->cd(2);
     histo_emu_E.Draw("");
   emuCan->cd(3);
     histo_emu_pt.Draw("");
  //Canvas for eTau
  TCanvas* etauCan = new TCanvas("etauCan", "Electron-Tau Events");
   etauCan->Divide(1,3);
   etauCan->cd(1);
     histo_etau_m.Draw("");
   etauCan->cd(2);
     histo_etau_E.Draw("");
   etauCan->cd(3);
     histo_etau_pt.Draw("");
  //Canvas for MuMu
  TCanvas* mumuCan = new TCanvas("mumuCan", "Muon-Muon Events");
   mumuCan->Divide(1,3);
   mumuCan->cd(1);
     histo_mumu_m.Draw("");
   mumuCan->cd(2);
     histo_mumu_E.Draw("");
   mumuCan->cd(3); 
     histo_mumu_pt.Draw("");
  //Canvas for MuTau
  TCanvas* mutauCan = new TCanvas("mutauCan", "Muon-Tau Events");
   mutauCan->Divide(1,3);
   mutauCan->cd(1);
     histo_mutau_m.Draw("");
   mutauCan->cd(2);
     histo_mutau_E.Draw("");
   mutauCan->cd(3);
     histo_mutau_pt.Draw("");
   //Canvas for TauTau
  TCanvas* tautauCan = new TCanvas("tautauCan", "Tau-Tau Events");
   tautauCan->Divide(1,3);
   tautauCan->cd(1);
     histo_tautau_m.Draw("");
   tautauCan->cd(2);
     histo_tautau_E.Draw("");
   tautauCan->cd(3);
     histo_tautau_pt.Draw("");

   //Canvas for specific events theta
  TCanvas* theta_spCan = new TCanvas("theta_spCan", "Theta Distro for Specific Events");
   theta_spCan->Divide(2,3);
   theta_spCan->cd(1);
     histo_ee_th.Draw("");
   theta_spCan->cd(2);
     histo_emu_th.Draw("");
   theta_spCan->cd(3);
     histo_etau_th.Draw("");
   theta_spCan->cd(4);
     histo_mumu_th.Draw("");
   theta_spCan->cd(5);
     histo_mutau_th.Draw("");
   theta_spCan->cd(6);
     histo_tautau_th.Draw("");


   //Canvas for specific events cos(theta)
  TCanvas* costh_spCan = new TCanvas("costh_spCan", "Cos(Theta) Distro for Specific Events");
   costh_spCan->Divide(2,3);
   costh_spCan->cd(1);
     histo_ee_costh.Draw("");
   costh_spCan->cd(2);
     histo_emu_costh.Draw("");
   costh_spCan->cd(3);
     histo_etau_costh.Draw("");
   costh_spCan->cd(4);
     histo_mumu_costh.Draw("");
   costh_spCan->cd(5);
     histo_mutau_costh.Draw("");
   costh_spCan->cd(6);
     histo_tautau_costh.Draw("");
  } 
  
 
  cout <<"\n\n nBins =  " <<nBins;
  cout << "\n nBin_sp = " <<nBins_sp ;
 
 
  cout << "\n\n Type Ctrl+c to return to prompt\n";
  graphicsPlease.Run();
 
 
 }


