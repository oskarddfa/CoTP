#include <TROOT.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraph.h>


void ROOTplot(int dim, double Phisquared[], double (*Potential)(double, double), double xstart, double xend, double height = 20)
	{

	
	auto xarray = new Double_t [dim];
	auto Varray = new Double_t [dim];
	Double_t dx = (xend-xstart)/(dim+1.0);
	
	for (int i=1; i<dim+1; i++)
		{
		xarray[i-1]= xstart + i*dx;
		Varray[i-1]= Potential(xarray[i-1], height);
		}
		
   TCanvas *c1 = new TCanvas("c1","Phi**2",200,10,700,500);
   c1->SetGrid();
   
   auto mgr = new TMultiGraph();
   
   auto Wellenfgr = new TGraph(dim,xarray,Phisquared);
   Wellenfgr->SetLineColor(2);
   Wellenfgr->SetLineWidth(4);
   Wellenfgr->SetMarkerColor(2);
   Wellenfgr->SetMarkerStyle(21);
   
   
   auto Potentialgr = new TGraph(dim,xarray,Varray);
   Potentialgr->SetLineColor(3);
   Potentialgr->SetLineWidth(4);
   Potentialgr->SetMarkerColor(3);
   Potentialgr->SetMarkerStyle(21);
   
   
   mgr->Add(Wellenfgr);
   mgr->Add(Potentialgr);
   
   


   
   mgr->Draw("APL");
   
   // TCanvas::Update() draws the frame, after which one can change it
   c1->Update();
   //c1->GetFrame()->SetBorderSize(12);
   c1->Modified();
   gROOT->GetListOfCanvases()->Draw();
   
   mgr->SetTitle("Phi**2, x, y");
   //mgr->GetXaxis()->SetLimits(xstart,xend);
   mgr->SetMinimum(0.);
   mgr->SetMaximum(0.5);

		
	}
		
