function[]=make_table(theta,phi,MM,d,nf)
  % nf : figure number
  Mc=-1i*MM; Mc=(Mc+Mc');
  pas=all(eig(Mc)>1e-8);  % check if all eigenvalues are positive
  if pas~=1
    disp('This is not a passive medium.');
    disp('Please provide appropriate material parameters');
    return;
  else
  end;
  
  [rss1, rps1, rsp1, rpp1, tss1, tps1, tsp1, tpp1]=fresnel_film_top(theta,phi,MM,d);
  [rss2, rps2, rsp2, rpp2, tss2, tps2, tsp2, tpp2]=fresnel_film_top(theta,phi+pi,MM,d);
  [rss3, rps3, rsp3, rpp3, tss3, tps3, tsp3, tpp3]=fresnel_film_bottom(theta,phi+pi,MM,d);
  
  f=figure(nf);
  uit = uitable(f);

  dt={'rss', rss1, rss2, rss3; 'rps', rps1, rps2, rps3;
     'rsp', rsp1, rsp2, rsp3; 'rpp', rpp1, rpp2, rpp3;
     'tss', tss1, tss2, tss3; 'tps', tps1, tps2, tps3;
     'tsp', tsp1, tsp2, tsp3; 'tpp', tpp1, tpp2, tpp3;};
  uit.Data=dt;
  uit.Position=[40 200 484 175];
  uit.ColumnName = {' ';'Incidence';'Reflection';'Transmission';};
  uit.RowName=[]; 
  uit.ColumnWidth={30,150,150,150};

  dim=[0.1 0.25 0.8 0.18];
  str1='For deriving SKL-1 and SKL-2, A) compare relations between rss, rps, rsp, rpp in incidence and reflection directions, and B) compare relations between tss, tps, tsp, tpp in incidence and transmission directions';
  annotation('textbox',dim,'String',str1);

  dim2=[0.1 0.1 0.8 0.1];
  str2='For deriving SKL-3, compare relations between coefficients in reflection and transmission directions';
  annotation('textbox',dim2,'String',str2);

  set(gca,'XColor', 'none','YColor','none');
  str3=strcat('(\theta,\phi)/\pi= (',num2str(theta/pi),' , ',num2str(phi/pi),')','    thickness, d\omega/c=',num2str(d));
  title(str3);
  
  %{
  %the following code works only for R2019a and later versions
  s1 = uistyle('BackgroundColor','yellow');
  s2 = uistyle('BackgroundColor','cyan');
  addStyle(uit,s1,'cell',[1:4,2:3]);
  addStyle(uit,s2,'cell',[5:8,3:4]);
  %}
  
  return;
