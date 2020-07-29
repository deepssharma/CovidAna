#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <math.h>
#include "TDatime.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

Long64_t nentries;
Long64_t nbytes;

void make_covid_tree(int, int, const char* , const char*, const char*);

void define_tree();

void read_input_files();

void clear_variables();

void init_tree();

void init_histos();

void ana_tree(int, int);

void fit_histos(int, int);

void draw_fits(const char*, const char*, const char*);

void write_output();

void return_month_year(const string&);

void get_means();

void set_timebins_array();

vector<double > hosp_cum[59];

char str[255];
char str2[255];
string str1;

string name; string name1; string title;
string Condition;
//char Condition[100];//sprintf(Condition,"Confirmed>-999 && State==\"%s\"",states[i]);

char delimiter=',';
string state;
string  junk;
string date;
string confirmed;
string deaths;
string recovered;
string active;
string incRate;
string tested;
string hospitalized;
string mortRate;
string testRate;
string hospRate;

char State[255];
int Confirmed=-999;
int Deaths=-999;
int Recovered=-999;
int Active=-999;
float IncRate=-999;
int Tested=-999;
int Hospitalized=-999;
float MortRate=-999;
float TestRate=-999;
float HospRate=-999;
TDatime *datetime;
int Days =-999;

int countlines=0;
int countfiles=0;
const int timebins =109;
const int meanbins =109;

double dailyD[59][10000];
double dailyC[59][10000];
double dailyH[59][10000];

double dailyH_C[59][10000];
double dailyD_H[59][10000];

double dailyH4_C[59][10000];
double dailyD21_H[59][10000];


float DD[59][10000];
int CC[59][10000];
int HH[59][10000];

double meanH1[59];
double meanH2[59];
double meanH3[59];

const char *plotstate[3] ={"Massachusetts","Georgia","Florida"};
//const char *plotstate[3] ={"Ohio","New York","Virginia"};
//char *plotstate[3];// ={"Alabama","Arizona","Virginia"};

int COLORS[9]  = {kBlack, kRed,kOrange,kSpring+5,kGreen,kCyan,kBlue,kMagenta+0,kPink+9};

double daysC[timebins];

TLegend *leg1;
TLegend *leg2;

TTree *MyTree;
TFile *fout;
TCanvas *c2;
TCanvas *c3;
TCanvas *c4;

TFile *fin;
TChain *T;

TH2D *hconfirmed[59];
TH2D *hdeaths[59];
TH2D *hhospital[59];
TH2D *hmortR[59];
TH2D *hhosptR[59];
TH2D *hCH;

TProfile *hC[59];
TProfile *hD[59];
TProfile *hH[59];
TProfile *hMR[59];
TProfile *hHR[59];

TGraphErrors *gr_deriv[59];
TGraphErrors *gr_deaths[59];
TGraphErrors *gr_dailyHC[59];
TGraphErrors *gr_dailyDH[59];

TH1D *hH_daily[59];
TH1D *hD_daily[59];
TH1D *hC_daily[59];
TH1D *hH_C[59];
TH1D *hD_H[59];
TH1D *hH4_C[59];
TH1D *hD21_H[59];
TH2D *hdailyHosp;

TF1 *fHC[59];
TF1 *fDH[59];
TF1 *fHC4[59];
TF1 *fDH21[59];


const char *states[59]={"Alabama","Alaska","Arizona","Arkansas","California","Colorado","Connecticut","Delaware","Diamond Princess","District of Columbia","Florida","Georgia","Grand Princess","Guam","Hawaii","Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana","Nebraska","Nevada","New Hampshire","New Jersey","New Mexico","New York","North Carolina","North Dakota","Ohio","Oklahoma","Oregon","Pennsylvania","Puerto Rico","Rhode Island","South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont","Virginia","Washington","West Virginia","Wisconsin","Wyoming","American Samoa","Northern Mariana Islands","Recovered","Virgin Islands"};

const char *states1[59]={"Alabama","Alaska","Arizona","Arkansas","California","Colorado","Connecticut","Delaware","DiamondPrincess","DistrictofColumbia","Florida","Georgia","GrandPrincess","Guam","Hawaii","Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana","Nebraska","Nevada","NewHampshire","NewJersey","NewMexico","NewYork","NorthCarolina","NorthDakota","Ohio","Oklahoma","Oregon","Pennsylvania","PuertoRico","RhodeIsland","SouthCarolina","SouthDakota","Tennessee","Texas","Utah","Vermont","Virginia","Washington","WestVirginia","Wisconsin","Wyoming","AmericanSamoa","NorthernMarianaIslands","Recovered","VirginIslands"};



void make_covid_tree(int fetch_data =0, int timelagH =4, int timelagD = 21, const char *st1="Massachusetts", const char *st2 ="Georgia",const char* st3 ="Florida")
{
  
  //sprintf(plotstate[0],"%s",st1);
  //sprintf(plotstate[1],"%s",st2);
  //sprintf(plotstate[2],"%s",st3);
  gErrorIgnoreLevel = kFatal;
  if(fetch_data==1)
    {
      cout<<"Reading the newly fetched data"<<endl;
      define_tree();
      
      read_input_files();

      cout<<" Total files Read are: " <<countfiles<<endl; 
      fout = new TFile("covidtree.root","RECREATE");
      fout->cd();
      MyTree->Write();
      fout->Close();
    }
  
  else if(fetch_data==0)
    {
      clear_variables();

      set_timebins_array();
      
      cout<<"Analyzing the already generated file"<<endl;
      T     = new TChain("MyTree");

      T->Add("./covidtree.root");      
      cout<<T->GetEntries()<<endl;
      //timebins = T->GetEntries/59;

      init_tree();
      nentries = T->GetEntries();
      nbytes = 0;
      
      init_histos();
      
      ana_tree(timelagH, timelagD);

      get_means();
      
      fit_histos(timelagH, timelagD);

      //cout<<st1<<"  "<<st2<<"  "<<st3<<endl;
      draw_fits(st1, st2, st3);
      
      write_output();
      
    }
}

void read_input_files()
{
  string inputlist=Form("data.list");
  //ifstream inputlist;
  //inputlist.open("data_full.list",ios::in);
  if(!inputlist.empty())//list of files
    {
      cout<<"Opening the list of files:"<<inputlist.c_str()<<endl;
      ifstream inputfilelist(inputlist.c_str());
      if(inputfilelist)//loop over file names
	{
	  while(inputfilelist)
	    {
	      inputfilelist.getline(str,255);
	      if(str[0]==0) break;
	      cout<<"Opening file: "<<str<<endl;
	      countfiles++;
	      //to get the date from the input file 
	      string dateN = string(str);
	      date =dateN.substr(84,85);
	      date.erase(10);
	      
	      return_month_year(date);
	      //cout<<datetime->GetDate()<<endl;
	      
	      if(str[0]!=0){
		ifstream filename(str);
		getline(filename, str1);
		//cout<<str1<<endl;
		countlines=0;
		while(filename)//loop over individual file
		  {
		    clear_variables();
		    //filename.getline(str2,255);
		    //filename>>state>>junk>>junk>>junk>>junk>>confirmed>>deaths>>recovered>>active>>junk>>incRate>>tested>>hospitalized>>mortRate>>junk>>junk>>testRate>>hospRate;
		    countlines++;
		    if(countlines==60) break;
		    getline(filename, state, delimiter);
		    getline(filename, junk, delimiter);
		    getline(filename, junk, delimiter);
		    getline(filename, junk, delimiter);
		    getline(filename, junk, delimiter);
		    getline(filename, confirmed, delimiter);
		    getline(filename, deaths, delimiter);
		    getline(filename, recovered, delimiter);
		    getline(filename, active, delimiter);
		    getline(filename, junk, delimiter);
		    getline(filename, incRate, delimiter);
		    getline(filename, tested, delimiter);
		    getline(filename, hospitalized, delimiter);
		    getline(filename, mortRate, delimiter);
		    getline(filename, junk, delimiter);
		    getline(filename, junk, delimiter);
		    getline(filename, testRate, delimiter);
		    getline(filename, hospRate, '\n');
		    
		    strcpy(State, state.c_str());
		    if(confirmed[0]!=0)Confirmed = stoi(confirmed);
		    if(deaths[0]!=0)Deaths    = stoi(deaths);
		    if(recovered[0]!=0)Recovered = stoi(recovered);
		    if(active[0]!=0)Active    = stoi(active);
		    if(incRate[0]!=0) IncRate  = std::stof(incRate);
		    if(tested[0]!=0)Tested    = stoi(tested);
		    if(hospitalized[0]!=0)Hospitalized=stoi(hospitalized);
		    if(mortRate[0]!=0) MortRate  = std::stof(mortRate);
		    if(testRate[0]!=0) TestRate  = std::stof(testRate);
		    if(hospRate.size()>1)HospRate  = std::stof(hospRate);
		    
		    MyTree->Fill();
		    //cout<<State<<"  "<<Recovered<<"  "<<Active<<"  "<<mortRate<<"  "<<HospRate<<endl;
		  }
	      }
	    }
	}//loop over file names
    }//list of files
}

void define_tree()
{
  MyTree = new TTree("MyTree","Covid-19 Tree");
  MyTree->Branch("State",&State,"State/C");
  MyTree->Branch("Confirmed",&Confirmed,"Confirmed/I");
  MyTree->Branch("Deaths",&Deaths,"Deaths/I");
  MyTree->Branch("Recovered",&Recovered,"Recovered/I");
  MyTree->Branch("Active",&Active,"Active/I");
  MyTree->Branch("IncRate",&IncRate,"IncRate/F");
  MyTree->Branch("Tested",&Tested,"Tested/I");
  MyTree->Branch("Hospitalized",&Hospitalized,"Hospitalized/I");
  MyTree->Branch("MortRate",&MortRate,"MortRate/F");
  MyTree->Branch("TestRate",&TestRate,"TestRate/F");
  MyTree->Branch("HospRate",&HospRate,"HospRate/F");
  MyTree->Branch("Days",&Days,"Days/I");
  MyTree->Branch("datetime",&datetime);
  
}

void return_month_year(const string& timestamp)// int day, int month, int year)
{
  int day, month, year, junk;
  year   = stoi(timestamp.substr(6));
  month  =stoi(timestamp.substr(0));
  day    =stoi(timestamp.substr(3));
  
  datetime->Set(year,month,day,0,0,0);
  TDatime dbegin;
  dbegin.Set(2020,04,12,0,0,0);
  int d1  = dbegin.Convert();
  int d2    =  datetime->Convert();
  Days  = d2/3600/24 - d1/3600/24;
  cout<<datetime->GetDate()<<"  "<<Days<<endl;
}

void clear_variables()
{
  memset(State,0,255);
  Confirmed=-999;
  Deaths=-999;
  Recovered=-999;
  Active=-999;
  IncRate=-999;
  Tested=-999;
  Hospitalized=-999;
  MortRate=-999;
  TestRate=-999;
  HospRate=-999;
  //Days =-999;

  for(int i =0; i<59; i++)
    {
      
      meanH1[i] =0;
      meanH2[i] =0;
      meanH3[i] =0;
      for(int j =0; j<10000; j++)
	{
	  dailyC[i][j]= 0;
	  dailyH[i][j]= 0;
	  dailyD[i][j]= 0;

	  dailyH_C[i][j]= 0;
	  dailyD_H[i][j]= 0;
	  dailyH4_C[i][j]= 0;
	  dailyD21_H[i][j]= 0;
	  
	  CC[i][j]  =0;
	  HH[i][j]  =0;
	  DD[i][j]  =0;

	  
	}
    }
}

void init_tree()
{
  T->SetBranchAddress("State",&State);
  T->SetBranchAddress("Confirmed",&Confirmed);
  T->SetBranchAddress("Deaths",&Deaths);
  T->SetBranchAddress("Recovered",&Recovered);
  T->SetBranchAddress("Active",&Active);
  T->SetBranchAddress("IncRate",&IncRate);
  T->SetBranchAddress("Tested",&Tested);
  T->SetBranchAddress("Hospitalized",&Hospitalized);
  T->SetBranchAddress("MortRate",&MortRate);
  T->SetBranchAddress("TestRate",&TestRate);
  T->SetBranchAddress("HospRate",&HospRate);
  T->SetBranchAddress("Days",&Days);
  T->SetBranchAddress("datetime",&datetime);

  
}

void ana_tree(int timelagH, int timelagD)
{
  double cum_number=0;
  for (int i=0; i<59; i++)
    {
      hC[i] = (TProfile*)hconfirmed[i]->ProfileX();
      hD[i] = (TProfile*)hdeaths[i]->ProfileX();
      hH[i] = (TProfile*)hhospital[i]->ProfileX();
      hMR[i]=(TProfile*)hmortR[i]->ProfileX();
      hHR[i]=(TProfile*)hhosptR[i]->ProfileX();
      //cum_num = hC[i]->GetBinContent(20);
      //hosp_cum.push_back
      //cout<<hC[i]->GetBinContent(48)<<endl;
    }
  
  int index=0; int st_index=-999;int cc=0;
  for (int i=0; i<nentries; i++)
  {
    
    nbytes += T->GetEntry(i);
    st_index=-999;
    cc = Days;
    if(i%59==0) { index=0;}
    else index++;
    //if(string(State)=="Virgin Islands") continue;
    //if(string(State)=="Massachusetts")cout<<Hospitalized<<"  "<<dailyD[st_index][cc]*100<<endl;
    //if(string(State)=="Georgia")cout<<Confirmed<<endl;
    for (int istate=0; istate<59; istate++)
      {
	if(string(State)==string(states[istate])) st_index=istate;
      }
    if(i<59)
      {
	CC[st_index][cc] = Confirmed;
	HH[st_index][cc] = Hospitalized;
	DD[st_index][cc] = Deaths;
	
      }
    else if(i>=59 && st_index!=-999){
      CC[st_index][cc] = Confirmed;
      HH[st_index][cc] = Hospitalized;
      DD[st_index][cc] = Deaths; 
            
      
      dailyC[st_index][cc] = Confirmed - CC[st_index][cc-1];
      dailyH[st_index][cc] = Hospitalized - HH[st_index][cc-1];
      dailyD[st_index][cc] = Deaths - DD[st_index][cc-1];
            
      if(dailyH[st_index][cc]!=0 && dailyC[st_index][cc]!=0)
	{
	  dailyH_C[st_index][cc] =  dailyH[st_index][cc]/dailyC[st_index][cc];
	  if(cc>=4 && dailyC[st_index][cc-timelagH]>0)dailyH4_C[st_index][cc] =  dailyH[st_index][cc]/dailyC[st_index][cc-timelagH];
	}
      if(dailyH[st_index][cc]!=0 && dailyD[st_index][cc]!=0)
	{
	  dailyD_H[st_index][cc] =  dailyD[st_index][cc]/dailyH[st_index][cc];
	  if(cc>=21 && dailyH[st_index][cc-timelagD]!=0)dailyD21_H[st_index][cc] =  dailyD[st_index][cc]/dailyH[st_index][cc-timelagD];
	}

      hhospital[st_index]->Fill(cc,dailyH[st_index][cc]);
    }
    //if(string(State)=="Massachusetts")cout<<datetime->GetDate()<<" "<<dailyD[st_index][cc]<<"  "<<dailyH[st_index][cc-14]<<"  "<<dailyD[st_index][cc-14]<<endl;
     
    //if(string(State)=="Georgia")cout<<dailyH[st_index][cc]<<",  "<<dailyD[st_index][cc]<<"  "<<dailyD_H[st_index][cc]<<endl;;
    }
  hdailyHosp=(TH2D*)hhospital[0]->Clone("hdailyHosp");
  for(int istate=0; istate<59; istate++)
    {
      /*gr_deriv[istate] = new TGraphErrors(meanbins, daysC, dailyH[istate],0,0);
	gr_deriv[istate]->SetMarkerColor(istate+1);
	gr_deriv[istate]->SetMarkerStyle(20);
	
	gr_deaths[istate] = new TGraphErrors(meanbins, daysC, dailyD[istate],0,0);
	gr_deaths[istate]->SetMarkerColor(istate+1);
	gr_deaths[istate]->SetMarkerStyle(20);
	
	gr_dailyHC[istate] = new TGraphErrors(timebins, daysC, dailyH_C[istate],0,0);
	gr_dailyHC[istate]->SetMarkerColor(istate+1);
	gr_dailyHC[istate]->SetMarkerStyle(24);
	
	gr_dailyDH[istate] = new TGraphErrors(timebins, daysC, dailyD_H[istate],0,0);
	gr_dailyDH[istate]->SetMarkerColor(istate+1);
	gr_dailyDH[istate]->SetMarkerStyle(20);
	
      */
      if(istate>0) 
	{
	  hdailyHosp->Add(hhospital[istate]);
	}
      for(int ibin=0; ibin<timebins; ibin++)
	  {
	    
	    hC_daily[istate]->SetBinContent(ibin, dailyC[istate][ibin]);
	    hH_daily[istate]->SetBinContent(ibin, dailyH[istate][ibin]);
	    hD_daily[istate]->SetBinContent(ibin, dailyD[istate][ibin]);
	    
	    hH_C[istate]->SetBinContent(ibin, dailyH_C[istate][ibin]);
	    hD_H[istate]->SetBinContent(ibin, dailyD_H[istate][ibin]);
	    
	    hH4_C[istate]->SetBinContent(ibin, dailyH4_C[istate][ibin]);
	    hD21_H[istate]->SetBinContent(ibin, dailyD21_H[istate][ibin]);
	    
	    //
	  }
      
	  
	  
    }
  //cout<<State<<"  "<<Days<<"  "<<dailyC[st_index][i]<<"  "<<dailyH[st_index][i]<<"  "<<dailyD[st_index][i]<<endl;
  
  //cout<<endl;
}

void init_histos()
{
 
  TH1F *htemp; TH1F *htemp1; TH1F *htemp2;
  int bin_hi =100; int nbins=100;
  
  for (int i=0; i<59; i++)
    {
      name = Form("Confirmed (%s)", states[i]);
      name1 = Form("hC_%s", states1[i]);
      Condition = Form("Confirmed>-999 && State==\"%s\"", states[i]);

      //T->Draw("Confirmed>>htemp()", Condition.c_str());
      T->Draw("Hospitalized>>htemp()", Condition.c_str());
      TH1F *h = (TH1F*)gROOT->FindObject("htemp");
      nbins   = h->GetNbinsX();
      bin_hi  = h->GetBinLowEdge(nbins);
      hconfirmed[i]=new TH2D(name1.c_str(), name.c_str(),timebins, 0, timebins, 200,0,bin_hi);
      //name=Form("Confirmed:Days>>%s",name1.c_str());
      name=Form("Hospitalized:Days>>%s",name1.c_str());
      T->Draw(name.c_str(),Condition.c_str());

      name = Form("Deaths (%s)", states[i]);
      name1 = Form("hD_%s", states1[i]);
      Condition = Form("Deaths>-999 && State==\"%s\"", states[i]);
      T->Draw("Deaths>>htemp1()", Condition.c_str());
      TH1F *h1 = (TH1F*)gROOT->FindObject("htemp1");
      nbins   = h1->GetNbinsX();
      bin_hi  = h1->GetBinLowEdge(nbins);
      hdeaths[i] = new TH2D(name1.c_str(),name.c_str(),timebins, 0, timebins, 200,0,bin_hi);
      name=Form("Deaths:Days>>%s",name1.c_str());
      T->Draw(name.c_str(),Condition.c_str());
      
      
      name = Form("Hospitalized (%s); Time since 04/12", states[i]);
      name1 = Form("hH_%s", states1[i]);
      Condition = Form("Hospitalized>-999 && State==\"%s\"", states[i]);
      T->Draw("Hospitalized>>htemp2()", Condition.c_str());
      TH1F *h2 = (TH1F*)gROOT->FindObject("htemp2");
      nbins   = h2->GetNbinsX();
      bin_hi  = h2->GetBinLowEdge(nbins);
      //hhospital[i]=new TH2D(name1.c_str(),name.c_str(),timebins, 0, timebins, 200,0,bin_hi);
      hhospital[i]=new TH2D(name1.c_str(),name.c_str(),timebins, 0, timebins, 8000,0,8000);
      name=Form("Hospitalized:Days>>%s",name1.c_str());
      //T->Draw(name.c_str(),Condition.c_str());
      //

      name = Form("Mortality Rate (%s)", states[i]);
      name1 = Form("hMortR_%s", states1[i]);
      Condition = Form("MortRate>-999 && State==\"%s\"", states[i]);
      hmortR[i] = new TH2D(name1.c_str(),name.c_str(),timebins, 0, timebins, 40,0,20);
      name=Form("MortRate:Days>>%s",name1.c_str());
      T->Draw(name.c_str(),Condition.c_str());

      name = Form("Hospitalized Rate (%s)", states[i]);
      name1 = Form("hHosptR_%s", states1[i]);
      Condition = Form("HospRate>-999 && State==\"%s\"", states[i]);
      hhosptR[i] = new TH2D(name1.c_str(),name.c_str(),timebins, 0, timebins, 60,0,30);
      name=Form("HospRate:Days>>%s",name1.c_str());
      T->Draw(name.c_str(),Condition.c_str());
      

      name1 = Form("hdailyC_%s; Time since 04/12", states1[i]);
      hC_daily[i]  = new TH1D(name1.c_str(), name1.c_str(), timebins, 0, timebins);
      name1 = Form("hdailyH_%s; Time since 04/12", states1[i]);
      hH_daily[i]  = new TH1D(name1.c_str(), name1.c_str(), timebins, 0, timebins);
      name1 = Form("hdailyD_%s; Time since 04/12", states1[i]);
      hD_daily[i]  = new TH1D(name1.c_str(), name1.c_str(), timebins, 0, timebins);
      name1 = Form("hhosp_conf_%s; Time since 04/12", states1[i]);
      hH_C[i]  = new TH1D(name1.c_str(), name1.c_str(), timebins, 0, timebins);
      name1    = Form("hdeath_hosp_%s; Time since 04/12", states1[i]);
      hD_H[i]  = new TH1D(name1.c_str(), name1.c_str(), timebins, 0, timebins);

      name1 = Form("hhosp4_conf_%s; Time since 04/12", states1[i]);
      hH4_C[i]  = new TH1D(name1.c_str(), name1.c_str(), timebins, 0, timebins);
      name1    = Form("hdeath21_hosp_%s; Time since 04/12", states1[i]);
      hD21_H[i]  = new TH1D(name1.c_str(), name1.c_str(), timebins, 0, timebins);
      
      
    }
  //hCH   = new TH2D("hCH","hCH", timebins, 0, timebins, 5000, 0, 5000);
}

void write_output()
{
  TFile *fgraphs = new TFile("./graphs.root","RECREATE");
  fgraphs->cd();
  
  for(int istate=0; istate<59; istate++)
      {
	/*name = Form("grH_%s", states1[istate]);
	gr_deriv[istate]->Write(name.c_str(), TObject::kOverwrite);
	name = Form("grD_%s", states1[istate]);
	gr_deaths[istate]->Write(name.c_str(), TObject::kOverwrite);
	name = Form("grHC_%s", states1[istate]);
	gr_dailyHC[istate]->Write(name.c_str(), TObject::kOverwrite);
	name = Form("grDH_%s", states1[istate]);
	gr_dailyDH[istate]->Write(name.c_str(), TObject::kOverwrite);*/
	
	hC_daily[istate]->Write();
	hH_daily[istate]->Write();
	hD_daily[istate]->Write();
	hH_C[istate]->Write();
	hD_H[istate]->Write();
	hH4_C[istate]->Write();
	hD21_H[istate]->Write();
	
	fHC[istate]->Write();
	fDH[istate]->Write();
	fHC4[istate]->Write();
	fDH21[istate]->Write();

	hhospital[istate]->Write();
	
      }
  
  hdailyHosp->Write();
  fgraphs->Close();
  c2->SaveAs("./figs/Trends_daily_hosp_conf_deaths.pdf");
  c3->SaveAs("./figs/Trends_lagged_daily_hosp_conf_deaths.pdf");
  c4->SaveAs("./figs/daily_conf_hosp_deaths_barchart.pdf");
}

void fit_histos(int timelagH, int timelagD)
{
  //TCanvas *c1 = new TCanvas();
  TF1 *fexpo  = new TF1("fexpo","expo", 0, timebins);
  
  for(int istate=0; istate<59; istate++)
    {
      name = Form("funcHC_%s",states1[istate]);
      fHC[istate] = new TF1(name.c_str(),"expo", 0, timebins);
      
      name = Form("funcDH_%s",states1[istate]);
      fDH[istate] = new TF1(name.c_str(),"expo", 0, timebins);

      name = Form("funcHC4_%s",states1[istate]);
      fHC4[istate] = new TF1(name.c_str(),"expo", 0, timebins);

      name = Form("funcDH21_%s",states1[istate]);
      fDH21[istate] = new TF1(name.c_str(),"expo", 0, timebins);
      
      hH_C[istate]->SetMinimum(0);
      hH_C[istate]->Fit("fexpo","QL","R", 2, timebins-2);
      
      //gSystem->ProcessEvents();
      fHC[istate]->SetParameters(fexpo->GetParameter(0), fexpo->GetParameter(1));

      fexpo  = new TF1("fexpo","expo", 0, timebins);
      hD_H[istate]->SetMinimum(0);
      hD_H[istate]->Fit("fexpo","Q","R", 2, timebins-2);
      fDH[istate]->SetParameters(fexpo->GetParameter(0), fexpo->GetParameter(1));

      fexpo  = new TF1("fexpo","expo", 0, timebins);
      hH4_C[istate]->SetMinimum(0);
      hH4_C[istate]->Fit("fexpo","QL","R", timelagH, timebins-2);
      fHC4[istate]->SetParameters(fexpo->GetParameter(0), fexpo->GetParameter(1));

      fexpo  = new TF1("fexpo","expo", 0, timebins);
      hD21_H[istate]->SetMinimum(0);
      hD21_H[istate]->Fit("fexpo","Q","R", timelagD, timebins-2);
      fDH21[istate]->SetParameters(fexpo->GetParameter(0), fexpo->GetParameter(1));

    }
}


void draw_fits(const char *st1, const char *st2, const char* st3)
{
  c2  = new TCanvas("c2","Hosp. and Mort. Rates",1000, 400);
  c3  =new TCanvas("c3","Hosp. and Mort. Rates (lagged)",1000, 400);
  c4   = new TCanvas("c4","c4", 1000,700);
  c4  ->Divide(3,3);
  
  c2->Divide(3,1);
  
  c3->Divide(3,1);
  
  char tmp[100];
  TH1D *hdummy   = new TH1D("hFits","; Days since 04/12; ",timebins, 0, timebins);
  hdummy->GetYaxis()->SetRangeUser(0,1.4);
  hdummy->SetStats(0);
  c2 ->cd(1);
  hdummy->Draw();
  c2 ->cd(2);
  hdummy->Draw();
  c2->cd(3);
  hdummy->Draw();
  
  c3 ->cd(1);
  hdummy->Draw();
  c3->cd(2);
  hdummy->Draw();
  c3->cd(3);
  hdummy->Draw();
  
  
  
  for(int istate=0; istate<59; istate++)
    {
      //if(meanH1[istate]>80)
      title = Form("Daily_Confirmed_%s",states1[istate]);
      hC_daily[istate]->SetTitle(title.c_str());
      title = Form("Daily_Hospitalized_%s",states1[istate]);
      hH_daily[istate]->SetTitle(title.c_str());
      title = Form("Daily_Deaths_%s",states1[istate]);
      hD_daily[istate]->SetTitle(title.c_str());
      
      name  =Form("%s",states[istate]);
      leg1  = new TLegend(0.099, 0.75, 0.70, 0.90);
      leg1  ->SetHeader(name.c_str(), "C");
      leg1  ->AddEntry(fHC[istate], "(Daily Hosp./Daily Conf)","l");
      leg1  ->AddEntry(fDH[istate], "(Daily Death./Daily Hosp)","l");
      leg1  ->SetTextColor(kBlack);
      leg1  ->SetTextSize(0.041);

      //name  =Form("%s",states[istate]);
      leg2  = new TLegend(0.099, 0.75, 0.70, 0.90);
      leg2  ->SetHeader(name.c_str(), "C");
      leg2  ->AddEntry((TObject*)0, "Lagged curves", "C");
      leg2  ->AddEntry(fHC[istate], "(Daily Hosp./Daily Conf)","l");
      leg2  ->AddEntry(fDH[istate], "(Daily Death./Daily Hosp)","l");
      leg2  ->SetTextColor(kBlack);
      leg2  ->SetTextSize(0.041);
	   
      if(strcmp(states[istate],st1)==0)
	 {
	   c2->cd(1);
	   fHC[istate]->Draw("same");
	   fDH[istate]->Draw("same");
	   fHC[istate]->SetLineColor(1);
	   fDH[istate]->SetLineColor(1);
	   fDH[istate]->SetLineStyle(7);
	   leg1->Draw();

	   
	   c3->cd(1);
	   fHC4[istate]->Draw("same");
	   fDH21[istate]->Draw("same");
	   fHC4[istate]->SetLineColor(1);
	   fDH21[istate]->SetLineColor(1);
	   fDH21[istate]->SetLineStyle(7);
	   leg2->Draw();

	   
	   c4->cd(1);
	   hC_daily[istate]->Draw("bar");
	   hC_daily[istate]->SetFillStyle(0);
	   hC_daily[istate]->SetLineColor(1);
	   
	   c4->cd(4);
	   hH_daily[istate]->Draw("bar");
	   hH_daily[istate]->SetFillStyle(0);
	   hH_daily[istate]->SetLineColor(1);
	   c4->cd(7);
	   hD_daily[istate]->Draw("bar");
	   hD_daily[istate]->SetFillStyle(0);
	   hD_daily[istate]->SetLineColor(1);
	   
	 }
      if(strcmp(states[istate],st2)==0)
      	{
	  c2->cd(2);
	  fHC[istate]->Draw("same");
	  fDH[istate]->Draw("same");
	  fHC[istate]->SetLineColor(2);
	  fDH[istate]->SetLineColor(2);
	  fDH[istate]->SetLineStyle(7);
	  leg1  ->SetTextColor(kRed);
	  leg1->Draw();
	  

	  c3->cd(2);
	  fHC4[istate]->Draw("same");
	  fDH21[istate]->Draw("same");
	  fHC4[istate]->SetLineColor(2);
	  fDH21[istate]->SetLineColor(2);
	  fDH21[istate]->SetLineStyle(7);
	  leg2  ->SetTextColor(kRed);
	  leg2->Draw();


	  c4->cd(2);
	  hC_daily[istate]->Draw("bar");
	  hC_daily[istate]->SetFillStyle(0);
	  hC_daily[istate]->SetLineColor(2);
	  c4->cd(5);
	  hH_daily[istate]->Draw("bar");
	  hH_daily[istate]->SetFillStyle(0);
	  hH_daily[istate]->SetLineColor(2);
	  c4->cd(8);
	  hD_daily[istate]->Draw("bar");
	  hD_daily[istate]->SetFillStyle(0);
	  hD_daily[istate]->SetLineColor(2);
	}
      if(strcmp(states[istate],st3)==0)
      	{
	  c2->cd(3);
	  fHC[istate]->Draw("same");
	  fDH[istate]->Draw("same");
	  fHC[istate]->SetLineColor(4);
	  fDH[istate]->SetLineColor(4);
	  fDH[istate]->SetLineStyle(7);
	  leg1  ->SetTextColor(kBlue);
	  leg1->Draw();
	  
	  c3->cd(3);
	  fHC4[istate]->Draw("same");
	  fDH21[istate]->Draw("same");
	  fHC4[istate]->SetLineColor(4);
	  fDH21[istate]->SetLineColor(4);
	  fDH21[istate]->SetLineStyle(7);
	  leg2  ->SetTextColor(kBlue);
	  leg2->Draw();

	  
	  c4->cd(3);
	  hC_daily[istate]->Draw("bar");
	  hC_daily[istate]->SetFillStyle(0);
	  hC_daily[istate]->SetLineColor(4);
	  c4->cd(6);
	  hH_daily[istate]->Draw("bar");
	  hH_daily[istate]->SetFillStyle(0);
	  hH_daily[istate]->SetLineColor(4);
	  c4->cd(9);
	  hD_daily[istate]->Draw("bar");
	  hD_daily[istate]->SetFillStyle(0);
	  hD_daily[istate]->SetLineColor(4);
	}
      
      //if((fHC[istate]->Integral(0,30)/30)>0.34)cout<<states[istate]<<"  "<<fHC[istate]->Integral(0,30)/30<<endl;
    }
  

}

void get_means()
{
  cout<<" ###########################################"<<endl;
  cout<<" States that have daily hospitlizations above 100 per day between Apri 12th to May 10th: "<<endl;
  cout<<" ###########################################"<<endl;
  

  for(int istate=0; istate<59; istate++)
    {
      hhospital[istate]->GetXaxis()->SetRangeUser(0,30);
      meanH1[istate] = hhospital[istate]->GetMean(2);
      if(meanH1[istate]>100)cout<<states[istate]<<"  "<<meanH1[istate]<<endl;
    }

  cout<<" ###########################################"<<endl<<endl;
  cout<<" States that have daily hospitlizations above 100 per day between May 11th to June 11th : "<<endl;
  cout<<" ###########################################"<<endl;
  for(int istate=0; istate<59; istate++)
  {
      hhospital[istate]->GetXaxis()->SetRangeUser(31,60);
      meanH2[istate] = hhospital[istate]->GetMean(2);
      if(meanH2[istate]>100)cout<<states[istate]<<"  "<<meanH2[istate]<<endl;
    }
  
  
  cout<<" ###########################################"<<endl<<endl;
  cout<<" States that have daily hospitlizations above 100 per day after June 12th: "<<endl;
  cout<<" ###########################################"<<endl;
  for(int istate=0; istate<59; istate++)
  {
      hhospital[istate]->GetXaxis()->SetRangeUser(61,106);
      meanH3[istate] = hhospital[istate]->GetMean(2);
      if(meanH3[istate]>100)cout<<states[istate]<<"  "<<meanH3[istate]<<endl;
    }
  
}

void set_timebins_array()
{
  for (int itime=0; itime<timebins; itime++)
    {
      daysC[itime]  = itime; 
    }
}
