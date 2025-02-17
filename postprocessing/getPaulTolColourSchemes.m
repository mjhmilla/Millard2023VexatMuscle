%%
% SPDX-FileCopyrightText: 2024 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%%
function cs = getPaulTolColourSchemes(nameOfColourScheme)
%%
% Returns the colour-blind safe colour schemes that are described in
% Paul Tol's blog entry:
% https://personal.sron.nl/~pault/#sec:qualitative
%%

switch nameOfColourScheme
    case 'bright'
        cs.blue     = [68,119,170]./255;
        cs.cyan     = [102,204,238]./255;
        cs.green    = [34,136,51]./255;
        cs.yellow   = [204,186,68]./255;
        cs.red      = [238,102,119]./255;
        cs.purple   = [170,51,119]./255;
        cs.grey     = [187,187,187]./255;
        
    case 'highContrast'
        cs.white    = [255,255,255]./255;
        cs.yellow   = [221,170,51]./255;
        cs.red      = [187,85,102]./255;
        cs.blue     = [0,68,136]./255;
        cs.black    = [0,0,0]./255;        
     
    case 'vibrant'
        cs.blue     = [0,119,187]./255;
        cs.cyan     = [51,187,238]./255;
        cs.teal     = [0,153,136]./255;
        cs.orange   = [238,119,51]./255;
        cs.red      = [204,51,17]./255;        
        cs.magenta  = [238,51,119]./255; 
        cs.grey     = [187,187,187]./255;

    case 'muted'
        cs.indigo   = [51,34,136]./255;
        cs.cyan     = [136,206,238]./255;
        cs.teal     = [68,170,153]./255;
        cs.green    = [17,119,51]./255;
        cs.olive    = [153,153,51]./255;
        cs.sand     = [221,204,119]./255;
        cs.rose     = [204,102,119]./255;        
        cs.wine     = [136,34,85]./255; 
        cs.purple   = [170,68,153]./255;        
        cs.palegrey = [221,221,221]./255;

    case 'mediumContrast'
        cs.white        =[255,255,255]./255;
        cs.lightYellow  =[238,204,102]./255;
        cs.lightRed     =[238,153,170]./255;
        cs.lightBlue    =[102,153,204]./255;
        cs.darkYellow   =[153,119,0]./255;        
        cs.darkRed      =[153,68,85]./255 ;
        cs.darkBlue     =[0,68,136]./255;
        cs.black        =[0,0,0]./255;

    case 'paleDark'
        cs.paleBlue     = [187,204,238]./255;
        cs.paleCyan     = [204,238,255]./255;
        cs.paleGreen    = [204,221,170]./255;
        cs.paleYellow   = [238,238,187]./255;
        cs.paleRed      = [255,204,204]./255;
        cs.paleGrey     = [221,221,221]./255;   

        cs.darkBlue     = [34,34,85]./255;
        cs.darkCyan     = [34,85,85]./255;
        cs.darkGreen    = [34,85,34]./255;
        cs.darkYellow   = [102,102,51]./255;
        cs.darkRed      = [102,51,51]./255;
        cs.darkGrey     = [85,85,85]./255;                
    
    case 'light'
        cs.lightBlue    = [119,170,221]./255;
        cs.lightCyan    = [153,221,255]./255;
        cs.mint         = [68,187,153]./255;
        cs.pear         = [187,204,51]./255;
        cs.olive        = [170,170,0]./255;
        cs.lightYellow  = [238,221,136]./255;
        cs.orange       = [238,136,102]./255;
        cs.pink         = [255,170,187]./255;
        cs.paleGrey     = [221,221,221]./255;

end



