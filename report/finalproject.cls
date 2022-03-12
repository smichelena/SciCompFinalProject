\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{finalproject}[30/10/21 Class for final project]

\LoadClass[onecolumn]{IEEEtran}

\RequirePackage{scalerel,stackengine}

\RequirePackage{parskip}
\setlength{\parindent}{0pt}

\RequirePackage{bbm}
\RequirePackage{float}
\RequirePackage{graphicx}
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage[margin=2cm]{geometry}
\RequirePackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc}
\RequirePackage{fancyhdr}

\usepackage[dvipsnames]{xcolor}
\definecolor{faded}{gray}{0.4}

\RequirePackage[numbered,framed]{matlab-prettifier}
\lstset{
  style              = Matlab-editor,
  basicstyle         = \mlttfamily,
  escapechar         = ",
  mlshowsectionrules = true,
}

%\RequirePackage{subfig}
%\RequirePackage[utf8]{inputenc}
%\RequirePackage{bm}
%\RequirePackage{mathtools}
%\RequirePackage{wrapfig}
%\RequirePackage{bigfoot} % to allow verbatim in footnote

\RequirePackage{physics}
\RequirePackage[section]{minted}
\definecolor{mintedbackground}{rgb}{0.95,0.95,0.95}

\newmintedfile[pythoncode]{Python}{
bgcolor=mintedbackground,
fontfamily=tt,
linenos=true,
numberblanklines=true,
numbersep=5pt,
gobble=0,
frame=leftline,
framerule=0.4pt,
framesep=2mm,
funcnamehighlighting=true,
tabsize=4,
obeytabs=false,
mathescape=false
samepage=false, %with this setting you can force the list to appear on the same page
showspaces=false,
showtabs =false,
texcl=false,
}


% Number subsections with letters since that's what professors do
% This prevents having like "1.1 Exercise part c: Solve for.." 
% And just allows "1.c  Solve for.."
\renewcommand{\thesubsection}{\thesection.\alph{subsection}}
% Use subection for equations (so equations match each answer)
\numberwithin{equation}{subsection}

% Commands
\newcommand{\pth}[1]{\left(#1\right)}
\newcommand{\bracket}[1]{\left[#1\right]}
\newcommand{\mtr}[1]{\begin{bmatrix}#1\end{bmatrix}}
\newcommand{\squiggly}[1]{\left\{#1\right\}}
%\newcommand{\norm}[1]{\left\Vert#1\right\Vert}
\def\ddfrac#1#2{\displaystyle\frac{\displaystyle #1}{\displaystyle #2}}


\stackMath
\newcommand\reallywidehat[1]{%
\savestack{\tmpbox}{\stretchto{%
  \scaleto{%
    \scalerel*[\widthof{\ensuremath{#1}}]{\kern-.6pt\bigwedge\kern-.6pt}%
    {\rule[-\textheight/2]{1ex}{\textheight}}%WIDTH-LIMITED BIG WEDGE
  }{\textheight}% 
}{0.5ex}}%
\stackon[1pt]{#1}{\tmpbox}%
}
