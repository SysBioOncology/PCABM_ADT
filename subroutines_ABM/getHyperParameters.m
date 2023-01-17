function sysVars=getHyperParameters(expname,template)

sysVars = template;

switch char(expname) %parameters change according to case name
      
    case 'abl' %only resistant tumor cells
        sysVars.params.TUpprol = 0.0596; 
        sysVars.params.TUpmig = 0.1167;
        sysVars.params.TUpmax = 50;
        
    case 'R1881' %growth in R1881
        sysVars.params.TUpprol = 0.1144; 
        sysVars.params.TUpmig = 0.1167;
        sysVars.params.M1pkill = 0.1116;
        sysVars.params.M2pkill = 0.0223;
        
    case 'DMSO' %growth without R1881
        sysVars.params.TUpprol = 0.0389;                 
        sysVars.params.TUpmig = 0.1;
        sysVars.params.M1pkill = 0.005;
        sysVars.params.M2pkill = 0.0349;
        sysVars.params.M2TUadd = 0.0995;
          
    case 'Castration_Resistance' %growth without R1881 with spontaneous resistance ability of tumor cells upon proliferation
        sysVars.params.TUpprol = 0.0389;                 %set lower proliferation probability
        sysVars.params.TUpmig = 0.1;
        sysVars.params.M1pkill = 0.005;
        sysVars.params.M2pkill = 0.0349;
        sysVars.params.M2TUadd = 0.0995;
        sysVars.params.TUpres = 0.001;                    
        sysVars.params.TUpprolres = 0.0596;             
        sysVars.params.TUpmigres = 0.1167;                
        sysVars.params.TUpmaxres = 50;
        
 
end

end