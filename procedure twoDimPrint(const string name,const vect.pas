procedure twoDimPrint(name :string;const vector<vector<long double> > & v)
begin
    WriteLn(name);
    if length>0 then
    begin
    for i:=0 to length-1 do begin
        for (int j = 0; j < v[i].size(); ++j) begin
            WriteLn(i,' ',j,' ',v[i][j]);
        end
    end
    end
    else
        WriteLn('is NULL');
end
procedure LinePrint(name :string;TLineSeries & v)
begin
    WriteLn(name);
    if length>0 then
    begin

    for i:=0 to length-1 do begin
            WriteLn(i,' ',v[i].first,' ',v[i].second);
    end
    end
    else
        WriteLn('is NULL');
end
procedure oneDimPrint(name :string;const vector<long double> & v)
begin
    WriteLn(name);
    if length>0 then
    begin

    for i:=0 to length-1 do begin
            WriteLn(i,' ',v[i]);
    end
    end
    else
        WriteLn('is NULL');
end

procedure mobilitySpectrum::logger()
begin
    WriteLn('Los started');
    WriteLn('A1',A1);
    twoDimPrint('Am',Am);
    WriteLn('An',An);
    WriteLn('B1',B1);
    oneDimPrint('B_spektr',B_spektr);
    WriteLn('Bn',Bn);
    twoDimPrint('Cl',Cl);
    twoDimPrint('Cl_t',Cl_t);
    twoDimPrint('Cm',Cm);
    twoDimPrint('Cm_t',Cm_t);
    WriteLn('Coef1',Coef1);
    WriteLn('Coef2',Coef2);
    twoDimPrint('Cr',Cr);
    twoDimPrint('Cr_t',Cr_t);
    WriteLn('F_s',F_s);
    WriteLn('GridPoints',GridPoints);
    oneDimPrint('GxxExp',GxxExp);
    oneDimPrint('Gxx_MC',Gxx_MC);
    oneDimPrint('Gxx_sp',Gxx_sp);
    oneDimPrint('GxyExp',GxyExp);
    oneDimPrint('Gxy_MC',Gxy_MC);
    oneDimPrint('Gxy_sp',Gxy_sp);

    oneDimPrint('IntGxx',IntGxx);
    oneDimPrint('IntGxy',IntGxy);
    oneDimPrint('IntMagField',IntMagField);
    oneDimPrint('Lv',Lv);

    WriteLn('MSLeft',MSLeft);
    WriteLn('MSRight',MSRight);
    oneDimPrint('MagField_spektr',MagField_spektr);
    WriteLn('MaxPoints',MaxPoints);
    WriteLn('Min_Spectr',Min_Spectr);

    oneDimPrint('Mobility',Mobility);

    WriteLn('Mu_max',Mu_max);
    WriteLn('Mu_min',Mu_min);
    oneDimPrint('Mv',Mv);

    WriteLn('NumberOfPoints',NumberOfPoints);
    WriteLn('PointPerInt',PointPerInt);
    WriteLn('Power_spektr',Power_spektr);
    twoDimPrint('Qm',Qm);
    WriteLn('SizeData',SizeData);
    oneDimPrint('Spectr_e',Spectr_e);
    oneDimPrint('Spectr_p',Spectr_p);
    WriteLn('Ves1',Ves1);
    WriteLn('Ves2',Ves2);
    oneDimPrint('Vpr',Vpr);
    WriteLn('W',W);
    oneDimPrint('Xr',Xr);
    oneDimPrint('Xv',Xv);
    WriteLn('bulua',bulua);
    LinePrint('electronMobilitySpectrum',electronMobilitySpectrum);
    LinePrint('holeMobilitySpectrum',holeMobilitySpectrum);

end