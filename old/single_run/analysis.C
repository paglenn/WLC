void analysis() {
	double delta = 0.001; 
	double PI = acos(-1);
	TCanvas *c = new TCanvas("Analysis","1e9 events");
	c->Divide(2,2);
	TFile * DataFile = new TFile("results.root","READ");

	vector<double> dcosx,dcosy,rpx,rpy,tpx,tpy,rptpx,rptpy;
	
	//dcos->Draw();
	for(int i = 0; i < 100; i++) {
		double y  = dcos->GetBinContent(i+1); 
		if(y>0) dcosx.push_back(dcos->GetBinCenter(i+1));
		if(y > 0) dcosy.push_back(log(y)) ;
		double y = rp->GetBinContent(i+1); 
		if(y>0) rpx.push_back(rp->GetBinCenter(i+1));
		if(y > 0) rpy.push_back(log(y)) ;
		double y = tp->GetBinContent(i+1);	
		if(y>0) tpx.push_back(tp->GetBinCenter(i+1));
		if(y > 0) tpy.push_back(log(y)) ;
		double y = rptp->GetBinContent(i+1);
		if(y>0) rptpx.push_back(rptp->GetBinCenter(i+1));
		if(y > 0) rptpy.push_back(log(y)) ;
	}

	c->cd(1);
	TGraph *tpG = new TGraph(tpx.size(),&tpx[0],&tpy[0]); 

	//data-driven gaussian fit
	double tp_mu = tp->GetMean();
	double tp_var = pow(tp->GetStdDev(),2) ; 
	std::cout<<"tp var: " << tp_var + pow(tp_mu,2)<< std::endl; 
	string tpString = "(-(x-[0])^2)/(2*[1]) - 0.5*log(2*[2]*[1])";
	TF1 *tpFunc = new TF1("tpGauss",tpString.c_str(),tpx[0],tpx[tpx.size()-1]); 
	tpFunc->SetParameters(tp_mu,tp_var,PI);

	// plot 
	tpG->SetTitle("T_{#perp}");
	tpG->Draw("AC*");
	tpFunc->Draw("same");
	tpG->GetYaxis()->SetTitle("probability dist with Gaussian curve");
	//tpFunc->Draw();

	c->cd(2);
	TGraph* rpG = new TGraph(rpx.size(),&rpx[0],&rpy[0]);

	// gaussian 
	double rp_mu = rp->GetMean();
	double rp_var = pow(rp->GetStdDev(),2) ; 
	std::cout<<"rp var: " << rp_var +pow(rp_mu,2) << std::endl; 
	string rpString = "(-(x-[0])^2)/(2*[1]) - 0.5*log(2*[2]*[1])";
	TF1 *rpFunc = new TF1("rpGauss",rpString.c_str(),rpx[0],rpx[rpx.size()-1]); 
	rpFunc->SetParameters(rp_mu,rp_var,PI);

	//plot 
	rpG->SetTitle("R_{#perp}");
	rpG->Draw("AC*");
	rpFunc->Draw("same");
	rpG->GetYaxis()->SetTitle("probability dist w/ Gaussian curve");
	//rpFunc->Draw();


	c->cd(3);
	TGraph* rptpG = new TGraph(rptpx.size(),&rptpx[0],&rptpy[0]);
	
	//data-driven gaussian fit
	double rptp_mu = rptp->GetMean();
	double rptp_var = pow(rptp->GetStdDev(),2.) ; 
	std::cout<<"rptp var: " << rptp_mu << std::endl; 
	string rptpString = "-(x-[0])*(x-[0])/(2.*[1]) - 0.5*log(2*[2]*[1])";
	TF1 *rptpFunc = new TF1("rptpGauss",rptpString.c_str(),rptpx[0],rptpx[rptpx.size()-1]); 
	rptpFunc->SetParameters(rptp_mu,rptp_var,PI);
	
	//plot 
	rptpG->SetTitle("R_{#perp} . T_{#perp}");
	rptpG->Draw("AC*");
	rptpG->GetYaxis()->SetTitle("probability dist w/ Gaussian curve");
	rptpFunc->Draw("same");
	//rptpFunc->Draw();

	//TKDE * rptpKDE = new KDE(

	c->cd(4); 
	TGraph* dcosG = new TGraph(dcosx.size(),&dcosx[0],&dcosy[0]);
	
	// distribution function
	string cosFuncString = "log([0])+[0]*(x-1)-log(1-exp(-2*[0]))";
	TF1 * cosFunc = new TF1("cos_dist",cosFuncString.c_str(),dcosx[0],dcosx[dcosx.size()-1]); 
	cosFunc->SetParameter(0,1./delta);
	
	// plot 
	dcosG->SetTitle("cos #theta");
	dcosG->Draw("AC*");
	dcosG->GetYaxis()->SetTitle("probability dist w/ theoretical curve");
	cosFunc->Draw("same");
	//cosFunc->Draw();

	//c->Update();
	//c->Print("graphs.pdf");
	//c->Modified();

		
	//exit(1);
}


