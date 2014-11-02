{ 

	Double_t L = 0.1; 
	TGraph * mygraph = new TGraphErrors("f(N).txt","%lg\t%lg"); 
	//mygraph.Draw(); 
	TF1 * myfunc = new TF1 ( "parametrization", " [0] + [1]*x", 55, 95); 
	myfunc->SetParameters(-25,-1); 

	TFitResultPtr fit = mygraph->Fit(myfunc, "M Q R S") ; 
	Double_t slope = fit->Value(1); 
	Double_t slope_err = fit->ParError(1); 
	Double_t PI = acos(-1); 
	double Fb_th  = PI*PI / ( 4 * L); 
	double Fb_meas = -1*slope; 
	double dF = slope_err; 
	cout << "expected: " << Fb_th  << endl; 
	cout << "measured: " << Fb_meas << " +/- " << fabs( 100*dF/Fb_meas) << "% " <<endl;

	mygraph->Draw("A 3"); 

}
