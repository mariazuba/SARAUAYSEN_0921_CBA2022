#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <LBPA.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  biolpar.allocate(1,8,"biolpar");
  nages.allocate("nages");
  nrep.allocate("nrep");
  pond.allocate(1,nrep,"pond");
  nlength.allocate("nlength");
  len_bins.allocate(1,nlength,"len_bins");
  LF_data.allocate(1,nrep,1,nlength,"LF_data");
  L50prior.allocate("L50prior");
  slopeprior.allocate("slopeprior");
  Fcrprior.allocate("Fcrprior");
  Loprior.allocate("Loprior");
  s1prior.allocate("s1prior");
  s2prior.allocate("s2prior");
  cv4.allocate("cv4");
  cv99.allocate("cv99");
  cv100.allocate("cv100");
  cv1.allocate("cv1");
  cv2.allocate("cv2");
  cv3.allocate("cv3");
  f3.allocate("f3");
  f4.allocate("f4");
  f2.allocate("f2");
  f7.allocate("f7");
  f5.allocate("f5");
  f6.allocate("f6");
 logL50ini=log(L50prior);
 logslopeini=log(slopeprior);
 logFcrini=log(Fcrprior);
 logLoini=log(Loprior);
 logs1ini=log(s1prior+1E-5);
 logs2ini=log(s2prior);
  ratio.allocate("ratio");
  h.allocate("h");
  nm.allocate("nm");
}

void model_parameters::initializationfunction(void)
{
  log_Fcr.set_initial_value(logFcrini);
  log_L50.set_initial_value(logL50ini);
  log_rango.set_initial_value(logslopeini);
  log_alfa.set_initial_value(logs1ini);
  log_beta.set_initial_value(logs2ini);
  log_Lo.set_initial_value(logLoini);
  log_Ftar.set_initial_value(-0.5);
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_Fcr.allocate(f2,"log_Fcr");
  log_L50.allocate(f3,"log_L50");
  log_rango.allocate(f4,"log_rango");
  log_alfa.allocate(f5,"log_alfa");
  log_beta.allocate(f6,"log_beta");
  log_Lo.allocate(f7,"log_Lo");
  log_Ftar.allocate(5,"log_Ftar");
  N0.allocate(1,nages,"N0");
  #ifndef NO_AD_INITIALIZE
    N0.initialize();
  #endif
  Ntar.allocate(1,nages,"Ntar");
  #ifndef NO_AD_INITIALIZE
    Ntar.initialize();
  #endif
  N.allocate(1,nages,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  Sel_a.allocate(1,nages,"Sel_a");
  #ifndef NO_AD_INITIALIZE
    Sel_a.initialize();
  #endif
  Sel.allocate(1,nlength,"Sel");
  #ifndef NO_AD_INITIALIZE
    Sel.initialize();
  #endif
  F.allocate(1,nages,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Z.allocate(1,nages,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(1,nages,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  mu_edad.allocate(1,nages,"mu_edad");
  #ifndef NO_AD_INITIALIZE
    mu_edad.initialize();
  #endif
  sigma_edad.allocate(1,nages,"sigma_edad");
  #ifndef NO_AD_INITIALIZE
    sigma_edad.initialize();
  #endif
  wmed.allocate(1,nlength,"wmed");
  #ifndef NO_AD_INITIALIZE
    wmed.initialize();
  #endif
  msex.allocate(1,nlength,"msex");
  #ifndef NO_AD_INITIALIZE
    msex.initialize();
  #endif
  Ps.allocate(1,nages,"Ps");
  #ifndef NO_AD_INITIALIZE
    Ps.initialize();
  #endif
  pred_Ctot_a.allocate(1,nages,"pred_Ctot_a");
  #ifndef NO_AD_INITIALIZE
    pred_Ctot_a.initialize();
  #endif
  pred_Ctot.allocate(1,nlength,"pred_Ctot");
  #ifndef NO_AD_INITIALIZE
    pred_Ctot.initialize();
  #endif
  likeval.allocate(1,10,"likeval");
  #ifndef NO_AD_INITIALIZE
    likeval.initialize();
  #endif
  edades.allocate(1,nages,"edades");
  #ifndef NO_AD_INITIALIZE
    edades.initialize();
  #endif
  prop_obs.allocate(1,nrep,1,nlength,"prop_obs");
  #ifndef NO_AD_INITIALIZE
    prop_obs.initialize();
  #endif
  prop_pred.allocate(1,nlength,"prop_pred");
  #ifndef NO_AD_INITIALIZE
    prop_pred.initialize();
  #endif
  Linf.allocate("Linf");
  #ifndef NO_AD_INITIALIZE
  Linf.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  Lo.allocate("Lo");
  #ifndef NO_AD_INITIALIZE
  Lo.initialize();
  #endif
  M.allocate("M");
  #ifndef NO_AD_INITIALIZE
  M.initialize();
  #endif
  s.allocate("s");
  #ifndef NO_AD_INITIALIZE
  s.initialize();
  #endif
  SPR.allocate("SPR");
  #ifndef NO_AD_INITIALIZE
  SPR.initialize();
  #endif
  SPRtar.allocate("SPRtar");
  #ifndef NO_AD_INITIALIZE
  SPRtar.initialize();
  #endif
  Fref.allocate("Fref");
  #ifndef NO_AD_INITIALIZE
  Fref.initialize();
  #endif
  YPR.allocate("YPR");
  #ifndef NO_AD_INITIALIZE
  YPR.initialize();
  #endif
  BPR.allocate("BPR");
  #ifndef NO_AD_INITIALIZE
  BPR.initialize();
  #endif
  Prob_talla.allocate(1,nages,1,nlength,"Prob_talla");
  #ifndef NO_AD_INITIALIZE
    Prob_talla.initialize();
  #endif
  FrecL.allocate(1,nages,1,nlength,"FrecL");
  #ifndef NO_AD_INITIALIZE
    FrecL.initialize();
  #endif
  alfa.allocate("alfa");
  #ifndef NO_AD_INITIALIZE
  alfa.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  dts.allocate("dts");
  #ifndef NO_AD_INITIALIZE
  dts.initialize();
  #endif
  B0.allocate("B0");
  #ifndef NO_AD_INITIALIZE
  B0.initialize();
  #endif
  slope.allocate("slope");
  #ifndef NO_AD_INITIALIZE
  slope.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
 Linf=biolpar(1);
 k=biolpar(2);
 M=biolpar(3);
 dts=biolpar(8);
 Ps=0.0;
 
}

void model_parameters::userfunction(void)
{
  f =0.0;
 Prob_length2age();
 Pop_Dynamic();
 Log_likelihood();
}

void model_parameters::Prob_length2age(void)
{
  int i, j;
 edades.fill_seqadd(1,1);
 mu_edad(1)=exp(log_Lo);
 for (i=2;i<=nages;i++)
  {
  mu_edad(i)=Linf*(1-exp(-k))+exp(-k)*mu_edad(i-1);
  }
  sigma_edad=exp(log_alfa)+exp(log_beta)*mu_edad;
  Prob_talla = ALK( mu_edad, sigma_edad, len_bins);
  slope=biolpar(7)-biolpar(6);
  wmed=exp(biolpar(4))*pow(len_bins,biolpar(5));
  msex=1./(1+exp(-log(19)*(len_bins-biolpar(6))/(slope)));
}

dvar_matrix model_parameters::ALK(dvar_vector& mu, dvar_vector& sig, dvector& x)
{
	//RETURN_ARRAYS_INCREMENT();
	int i, j;
	dvariable z1;
	dvariable z2;
	int si,ni; si=mu.indexmin(); ni=mu.indexmax();
	int sj,nj; sj=x.indexmin(); nj=x.indexmax();
	dvar_matrix pdf(si,ni,sj,nj);
	pdf.initialize();
	double xs=0.5*(x[sj+1]-x[sj]);
	for(i=si;i<=ni;i++) //loop over ages
	{
		 for(j=sj;j<=nj;j++) //loop over length bins
		{
			z1=((x(j)-xs)-mu(i))/sig(i);
			z2=((x(j)+xs)-mu(i))/sig(i);
			pdf(i,j)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
		pdf(i)/=sum(pdf(i));
	}//end nage
	//RETURN_ARRAYS_DECREMENT();
	return(pdf);
}

void model_parameters::Pop_Dynamic(void)
{
  Sel=1./(1+exp(-log(19)*(len_bins-exp(log_L50))/(exp(log_rango))));
  Sel_a=Prob_talla*Sel;
  F=exp(log_Fcr)*Sel_a;
  Z=F+M;
  S=exp(-1.*Z);
 // Unfished biomass////////////////////////////
  N0(1)=1.0;
  for (int j=2;j<=nages;j++)
  { N0(j)=N0(j-1)*exp(-1.*M);}
    N0(nages)=N0(nages)/(1-exp(-1.*M));
  B0=sum(elem_prod((N0*exp(-dts*M))*Prob_talla,elem_prod(wmed,msex)));
  alfa=4*h/(5*h-1);
  beta=(1-h)/(5*h-1)*B0;
 // LF estimation ////////////////////////////
  N(1)=1.0;
  for (int i=2;i<=nages;i++){
  N(i)=N(i-1)*exp(-Z(i-1));
  }
  N(nages)=N(nages)/(1-exp(-Z(nages)));
  pred_Ctot_a=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
  pred_Ctot=pred_Ctot_a*Prob_talla;
 // Proportions ////////////////////////////
  for (int j=1;j<=nrep;j++){
  prop_obs(j)=LF_data(j)/sum(LF_data(j));
  }
  prop_pred=pred_Ctot/sum(pred_Ctot);
 // SPR estimation and Ftarg////////////////////////////
  SPR=1/B0*(alfa*sum(elem_prod(elem_prod(N,exp(-dts*Z))*Prob_talla,elem_prod(wmed,msex)))-beta);
  Ntar(1)=1.0;
  Z=M+exp(log_Ftar)*Sel_a;
  for (int i=2;i<=nages;i++){
  Ntar(i)=Ntar(i-1)*exp(-Z(i-1));
  }
  Ntar(nages)=Ntar(nages)/(1-exp(-Z(nages)));
  // Per-recruit biomass and yield /////////////////////////
  SPRtar=1/B0*(alfa*sum(elem_prod(elem_prod(Ntar,exp(-dts*Z))*Prob_talla,elem_prod(wmed,msex)))-beta);
}

void model_parameters::Log_likelihood(void)
{
  s=0;
  for (int j=1;j<=nrep;j++){
  s+=-nm*sum(pond(j)*elem_prod(prop_obs(j),log(prop_pred)));// LF_data
  }
  likeval(1)=s;// LF_data
  likeval(2)=0.5*square((log_Lo-logLoini)/cv1);
  likeval(3)=0.5*square((log_alfa-logs1ini)/cv2);
  likeval(4)=0.5*square((log_beta-logs2ini)/cv3);
  likeval(5)=0.5*square((log_L50-logL50ini)/cv4);
  likeval(6)=0.5*square((log_rango-logslopeini)/cv99);
  likeval(7)=0.5*square((log_Fcr-logFcrini)/cv100);
  if(last_phase()){
  likeval(8)=1000*square((log(SPRtar)-log(ratio)));}
  f=sum(likeval);
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  for (int i=1;i<=nages;i++){
    FrecL(i)=Prob_talla(i)*pred_Ctot_a(i);
  }
  report << "Length_bins" << endl;
  report <<  len_bins << endl;
  report << "Observed_frequencies" << endl;
  report <<  prop_obs << endl;
  report << "Predicted_frequency" << endl;
  report <<  prop_pred << endl;
  report << "Catch_length_frequency (columns) by age groups (rows)" << endl;
  report <<  FrecL  << endl;
  report << "Probability_of_length (columns) at-age (rows)" << endl;
  report <<  Prob_talla  << endl;
  report << "Age_Length_s.e_N_Catch_Selectivity_F_Weight_Maturity " << endl;
  for (int j=1;j<=nages;j++){ // l
  report << edades(j) <<" "<<mu_edad(j)<<" "<<sigma_edad(j)<<" "<<N(j)<<" "<<pred_Ctot_a(j)<<" "<<Sel_a(j)<<" "<<Sel_a(j)*exp(log_Fcr)<<" "<<Prob_talla(j)*wmed<<" "<<Prob_talla(j)*msex<<endl; 
  }
    for (int i=1;i<=nages;i++){
    if(mu_edad(i)>=biolpar(1)*2/3){
    Ps(i)=1;
    }};
  report << "Length_frequency_of_exploitable_population_current_target_unfished" << endl;
  report << elem_prod(N,Sel_a)*Prob_talla<<endl;
  report << elem_prod(Ntar,Sel_a)*Prob_talla<<endl;
  report << elem_prod(N0,Sel_a)*Prob_talla<<endl;
  report << "Selectivity_and_maturity_at_length" << endl;
  report << Sel<<endl;
  report << msex<<endl;
  report << "Model parameters " << endl;
  report<<"F_L50_slope_a0_cv_Lr_Ftar"<<endl;
  report<<exp(log_Fcr)<<" "<<exp(log_L50)<<" "<<exp(log_rango)<<" "<<exp(log_alfa)<<" "<<exp(log_beta)<<" "<<exp(log_Lo)<<" "<<exp(log_Ftar)<<endl;
  report<<"F/Ftar_SPR_SPRtar"<<endl;
  report<<exp(log_Fcr)/exp(log_Ftar)<<" "<<SPR<<" "<<SPRtar<<endl;
  report << "Log-likelihood_components" << endl;
  report << "Proportions" << endl;
  report << likeval(1) << endl;
  report << "Lr" << endl;
  report << likeval(2) << endl;
  report << "a0" << endl;
  report << likeval(3) << endl;
  report << "cv" << endl;
  report << likeval(4) << endl;
  report << "L50" << endl;
  report << likeval(5) << endl;
  report << "slope" << endl;
  report << likeval(6) << endl;
  report << "F" << endl;
  report << likeval(7) << endl;
  report << "Total" << endl;
  report << sum(likeval) << endl;
  report <<"Per_recruit_Analysis"<<endl;
  report <<"F_Y/R_SSB/R"<<endl;
   Fref=0.0;
   while(Fref<=3*exp(log_Fcr)){
    F=Fref*Sel_a;
    Z=F+M;
    S=exp(-1.*Z);
    N(1)=1.0;
      for (int i=2;i<=nages;i++){
    N(i)=N(i-1)*exp(-Z(i-1));
    }
  N(nages)=N(nages)/(1-exp(-Z(nages)));
  BPR=alfa*sum(elem_prod(elem_prod(N,exp(-dts*Z))*Prob_talla,elem_prod(wmed,msex)))-beta;
  YPR=(alfa*BPR/(beta+BPR))*sum(elem_prod(elem_prod(elem_div(F,Z),elem_prod(1.-S,N))*Prob_talla,wmed));
  report <<Fref<<" "<<YPR<<" "<<BPR<<endl;
  Fref+=0.01;
  }
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
