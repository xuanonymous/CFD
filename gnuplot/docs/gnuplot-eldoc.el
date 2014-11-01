;;;; Automatically generated by doc2texi.el on Fri, 08 March 2013

;;;; This file is generated automatically from the Gnuplot
;;;; documentation by `doc2texi.el', part of the Gnuplot distribution.
;;;; It is not intended to be edited manually.  See docs/gnuplot.doc
;;;; in the gnuplot source tree (available from
;;;; gnuplot.sourcefourge.net) for the original.

;; This file is covered by the Gnuplot licensing terms:

;; Copyright 1986 - 1993, 1998, 2004   Thomas Williams, Colin Kelley
;;
;; Permission to use, copy, and distribute this software and its
;; documentation for any purpose with or without fee is hereby granted,
;; provided that the above copyright notice appear in all copies and
;; that both that copyright notice and this permission notice appear
;; in supporting documentation.
;;
;; Permission to modify the software is granted, but not the right to
;; distribute the complete modified source code.  Modifications are to
;; be distributed as patches to the released version.  Permission to
;; distribute binaries produced by compiling modified sources is granted,
;; provided you
;;   1. distribute the corresponding source modifications from the
;;    released version in the form of a patch file along with the binaries,
;;   2. add special version identification to distinguish your version
;;    in addition to the base release version number,
;;   3. provide your name and address as the primary contact for the
;;    support of your modified version, and
;;   4. retain our contact information in regard to use of the base
;;    software.
;; Permission to distribute the released version of the source code along
;; with corresponding source modifications in the form of a patch file is
;; granted with same provisions 2 through 4 for binary distributions.
;;
;; This software is provided "as is" without express or implied warranty
;; to the extent permitted by applicable law.

(eval-when-compile (defvar gnuplot-eldoc-hash nil))
(setq gnuplot-eldoc-hash (let ((tbl (make-hash-table :test (quote equal))) (alist (quote (("x11" "set terminal x11 {<n> | window \"<string>\"} [more ...]" "set terminal x11 {<n> | window \"<string>\"}
                 {title \"<string>\"}
                 {{no}enhanced} {font <fontspec>}
                 {linewidth LW} {solid|dashed}
                 {{no}persist} {{no}raise} {{no}ctrlq}
                 {close}
                 {size XX,YY} {position XX,YY}
set terminal x11 {reset}") ("wxt" "set term wxt {<n>} [more ...]" "set term wxt {<n>}
             {size <width>,<height>} {background <rgb_color>}
             {{no}enhanced}
             {font <font>} {fontscale <scale>}
             {title \"title\"}
             {dashed|solid} {dashlength <dl>}
             {{no}persist}
             {{no}raise}
             {{no}ctrl}
             {close}") ("windows" "set terminal windows {<n>} [more ...]" "set terminal windows {<n>}
                     {color | monochrome}
                     {solid | dashed}
                     {enhanced | noenhanced}
                     {font <fontspec>}
                     {fontscale <scale>}
                     {linewdith <scale>}
                     {background <rgb color>}
                     {title \"Plot Window Title\"}
                     {size <width>,<height>}
                     {position <x>,<y>}
                     {close}") ("vgagl" "set terminal vgagl \\\\ [more ...]" "set terminal vgagl \\\\
             background [red] [[green] [blue]] \\\\
             [uniform | interpolate] \\\\
             [dump \"file\"] \\\\
             [mode]") ("tpic" "set terminal tpic <pointsize> <linewidth> <interval>") ("tgif" "set terminal tgif {portrait | landscape | default} {<[x,y]>} [more ...]" "set terminal tgif {portrait | landscape | default} {<[x,y]>}
                  {monochrome | color}
                  {{linewidth | lw} <LW>}
                  {solid | dashed}
                  {font \"<fontname>{,<fontsize>}\"}") ("svg" "set terminal svg {size <x>,<y> {|fixed|dynamic}} [more ...]" "set terminal svg {size <x>,<y> {|fixed|dynamic}}
                 {{no}enhanced}
                 {fname \"<font>\"} {fsize <fontsize>}
                 {mouse} {standalone | jsdir <dirname>}
                 {name <plotname>}
                 {font \"<fontname>{,<fontsize>}\"}
                 {fontfile <filename>}
                 {rounded|butt} {solid|dashed} {linewidth <lw>}
                 {background <rgb_color>}") ("regis" "set terminal regis {4 | 16}\"") ("pstricks" "set terminal pstricks {hacktext | nohacktext} {unit | nounit}") ("pdf" "set terminal pdf {monochrome|color|colour} [more ...]" "set terminal pdf {monochrome|color|colour}
                 {{no}enhanced}
                 {fname \"<font>\"} {fsize <fontsize>}
                 {font \"<fontname>{,<fontsize>}\"} {fontscale <scale>}
                 {linewidth <lw>} {rounded|butt}
                 {solid|dashed} {dl <dashlength>}}
                 {size <XX>{unit},<YY>{unit}}") ("pbm" "set terminal pbm {<fontsize>} {<mode>} {size <x>,<y>}") ("Openstep_(next)" "set terminal openstep {<mode>} {<type> } {<color>} {<dashed>} [more ...]" "set terminal openstep {<mode>} {<type> } {<color>} {<dashed>}
           {\"<fontname>\"} {<fontsize>} title {\"<newtitle>\"}") ("next" "set terminal next {<mode>} {<type> } {<color>} {<dashed>} [more ...]" "set terminal next {<mode>} {<type> } {<color>} {<dashed>}
           {\"<fontname>\"} {<fontsize>} title {\"<newtitle>\"}") ("mif" "set terminal mif {color | colour | monochrome} {polyline | vectors} [more ...]" "set terminal mif {color | colour | monochrome} {polyline | vectors}
                 {help | ?}") ("mp" "set term mp {color | colour | monochrome} [more ...]" "set term mp {color | colour | monochrome}
            {solid | dashed}
            {notex | tex | latex}
            {magnification <magsize>}
            {psnfss | psnfss-version7 | nopsnfss}
            {prologues <value>}
            {a4paper}
            {amstex}
            {\"<fontname> {,<fontsize>}\"} ") ("macintosh" "set terminal macintosh {singlewin | multiwin} {vertical | novertical} [more ...]" "set terminal macintosh {singlewin | multiwin} {vertical | novertical}
                       {size <width>, <height> | default}") ("lua" "set terminal lua <target name> | \"<file name>\" [more ...]" "set terminal lua <target name> | \"<file name>\"
                    {<script_args> ...}
                    {help}") ("latex" "set terminal {latex | emtex} {default | {courier|roman} {<fontsize>}} [more ...]" "set terminal {latex | emtex} {default | {courier|roman} {<fontsize>}}
             {size <XX>{unit}, <YY>{unit}} {rotate | norotate}") ("imagen" "set terminal imagen {<fontsize>} {portrait | landscape} [more ...]" "set terminal imagen {<fontsize>} {portrait | landscape}
                    {[<horiz>,<vert>]}") ("hppj" "set terminal hppj {FNT5X9 | FNT9X17 | FNT13X25}") ("hpljii" "set terminal hpljii | hpdj {<res>}") ("hpgl" "set terminal pcl5 {mode <mode>} {<plotsize>} [more ...]" "set terminal pcl5 {mode <mode>} {<plotsize>}
    {{color {<number_of_pens>}} | monochrome} {solid | dashed}
    {font <font>} {size <fontsize>} {pspoints | nopspoints}") ("hpgl" "set terminal hpgl {<number_of_pens>} {eject}") ("hp500c" "set terminal hp500c {<res>} {<comp>}") ("gpic" "set terminal gpic {<x> <y>}") ("ggi" "set terminal ggi [acceleration <integer>] [[mode] {mode}]") ("png_" "set terminal png  [more ...]" "set terminal png 
       {{no}enhanced}
       {{no}transparent} {{no}interlace}
       {{no}truecolor} {rounded|butt}
       {linewidth <lw>} {dashlength <dl>}
       {tiny | small | medium | large | giant}
       {font \"<face> {,<pointsize>}\"} {fontscale <scale>}
       {size <x>,<y>} {{no}crop}
       {background <rgb_color>}") ("fig" "set terminal fig {monochrome | color} [more ...]" "set terminal fig {monochrome | color}
                 {landscape | portrait}
                 {small | big | size <xsize> <ysize>}
                 {metric | inches}
                 {pointsmax <max_points>}
                 {solid | dashed}
                 {font \"<fontname>{,<fontsize>}\"}
                 {textnormal | {textspecial texthidden textrigid}}
                 {{thickness|linewidth} <units>}
                 {depth <layer>}
                 {version <number>}") ("epson_180dpi" "set terminal dpu414 {small | medium | large} {normal | draft}") ("epson_180dpi" "set terminal nec_cp6 {monochrome | colour | draft}") ("emxvga" "set terminal emxvga [more ...]" "set terminal emxvga
set terminal emxvesa {vesa-mode}
set terminal vgal") ("emf" "set terminal emf {color | monochrome} {solid | dashed} [more ...]" "set terminal emf {color | monochrome} {solid | dashed}
                 {enhanced {noproportional}}
                 {rounded | butt}
                 {linewidth <LW>} {dashlength <DL>}
                 {size XX,YY} {background <rgb_color>}
                 {font \"<fontname>{,<fontsize>}\"}
                 {fontscale <scale>}") ("eepic" "set terminal eepic {default} {color|dashed} {rotate} {size XX,YY} [more ...]" "set terminal eepic {default} {color|dashed} {rotate} {size XX,YY}
                   {small|tiny|<fontsize>}") ("dumb" "set terminal dumb {size <xchars>,<ychars>} {[no]feed} [more ...]" "set terminal dumb {size <xchars>,<ychars>} {[no]feed}
                  {[no]enhanced}") ("svga" "set terminal svga {\"<fontname>\"}\"") ("corel" "set terminal corel {  default [more ...]" "set terminal corel {  default
                    | {monochrome | color
                         {\"<font>\" {<fontsize> 
                            {<xsize> <ysize> {<linewidth> }}}}}") ("context" "set term context {default} [more ...]" "set term context {default}
        {defaultsize | size <scale> | size <xsize>{in|cm}, <ysize>{in|cm}}
        {input | standalone}
        {timestamp | notimestamp}
        {noheader | header \"<header>\"}
        {color | colour | monochrome}
        {rounded | mitered | beveled} {round | butt | squared}
        {dashed | solid} {dashlength | dl <dl>}
        {linewidth | lw <lw>}
        {fontscale <fontscale>}
        {mppoints | texpoints}
        {inlineimages | externalimages}
        {defaultfont | font \"{<fontname>}{,<fontsize>}\"}") ("cgm" "set terminal cgm {color | monochrome} {solid | dashed} {{no}rotate} [more ...]" "set terminal cgm {color | monochrome} {solid | dashed} {{no}rotate}
                 {<mode>} {width <plot_width>} {linewidth <line_width>}
                 {font \"<fontname>,<fontsize>\"}
                 {background <rgb_color>}") ("canvas" "set terminal canvas {size <xsize>, <ysize>} {background <rgb_color>} [more ...]" "set terminal canvas {size <xsize>, <ysize>} {background <rgb_color>}
                    {font {<fontname>}{,<fontsize>}} | {fsize <fontsize>}
                    {{no}enhanced} {linewidth <lw>}
                    {rounded | butt}
                    {solid | dashed {dashlength <dl>}}
                    {standalone {mousing} | name '<funcname>'}
                    {jsdir 'URL/for/javascripts'}
                    {title '<some string>'}") ("be" "set terminal be {reset} {<n>}") ("aqua" "set terminal aqua {<n>} {title \"<wintitle>\"} {size <x> <y>} [more ...]" "set terminal aqua {<n>} {title \"<wintitle>\"} {size <x> <y>}
                  {font \"<fontname>{,<fontsize>}\"}
                  {{no}enhanced} {solid|dashed} {dl <dashlength>}}") ("aifm" "set terminal aifm {color|monochrome} {\"<fontname>\"} {<fontsize>}") ("While" "while (<expr>) { [more ...]" "while (<expr>) {
    <commands>
}") ("update" "update <filename> {<filename>}") ("test" "test {terminal | palette [rgb|rbg|grb|gbr|brg|bgr]}") ("stats_(Statistical_Summary)" "stats 'filename' [using N[:M]] [name 'prefix'] [[no]output]]") ("data-file" "splot '<file_name>' {binary <binary list>} [more ...]" "splot '<file_name>' {binary <binary list>}
                    {{nonuniform} matrix}
                    {index <index list>}
                    {every <every list>}
                    {using <using list>}") ("splot" "splot {<ranges>} [more ...]" "splot {<ranges>}
      {<iteration>}
      <function> | \"<datafile>\" {datafile-modifiers}}
      {<title-spec>} {with <style>}
      {, {definitions{,}} <function> ...}") ("zeroaxis" "set {x|x2|y|y2|z}zeroaxis { {linestyle | ls <line_style>} [more ...]" "set {x|x2|y|y2|z}zeroaxis { {linestyle | ls <line_style>}
                           | { linetype | lt <line_type>}
                             { linewidth | lw <line_width>}}
unset {x|x2|y|y2|z}zeroaxis
show {x|y|z}zeroaxis") ("zero" "set zero <expression> [more ...]" "set zero <expression>
show zero") ("xyplane" "set xyplane at <zvalue> [more ...]" "set xyplane at <zvalue>
set xyplane relative <frac>
set ticslevel <frac>        # equivalent to set xyplane relative
show xyplane") ("xtics" "set xtics {axis | border} {{no}mirror} [more ...]" "set xtics {axis | border} {{no}mirror}
          {in | out} {scale {default | <major> {,<minor>}}}
          {{no}rotate {by <ang>}} {offset <offset> | nooffset}
          {left | right | center | autojustify}
          {add}
          {  autofreq
           | <incr>
           | <start>, <incr> {,<end>}
           | ({\"<label>\"} <pos> {<level>} {,{\"<label>\"}...) }
          { format \"formatstring\" } { font \"name{,<size>}\" }
          { rangelimited }
          { textcolor <colorspec> }
unset xtics
show xtics") ("xrange" "set xrange { [{{<min>}:{<max>}}] {{no}reverse} {{no}writeback} } [more ...]" "set xrange { [{{<min>}:{<max>}}] {{no}reverse} {{no}writeback} }
           | restore
show xrange") ("xmtics" "set xmtics [more ...]" "set xmtics
unset xmtics
show xmtics") ("xlabel" "set xlabel {\"<label>\"} {offset <offset>} {font \"<font>{,<size>}\"} [more ...]" "set xlabel {\"<label>\"} {offset <offset>} {font \"<font>{,<size>}\"}
           {textcolor <colorspec>} {{no}enhanced}
           {rotate by <degrees> | rotate parallel | norotate}
show xlabel") ("xdtics" "set xdtics [more ...]" "set xdtics
unset xdtics
show xdtics") ("xdata" "set xdata {time} [more ...]" "set xdata {time}
show xdata") ("view" "set view <rot_x>{,{<rot_z>}{,{<scale>}{,<scale_z>}}} [more ...]" "set view <rot_x>{,{<rot_z>}{,{<scale>}{,<scale_z>}}}
set view map
set view {no}equal {xy|xyz}
show view") ("version" "show version {long}") ("variables" "show variables      # show variables that do not begin with GPVAL_ [more ...]" "show variables      # show variables that do not begin with GPVAL_
show variables all  # show all variables including those beginning GPVAL_
show variables NAME # show only variables beginning with NAME") ("title_" "set title {\"<title-text>\"} {offset <offset>} {font \"<font>{,<size>}\"} [more ...]" "set title {\"<title-text>\"} {offset <offset>} {font \"<font>{,<size>}\"}
          {{textcolor | tc} {<colorspec> | default}} {{no}enhanced}
show title") ("timefmt" "set timefmt \"<format string>\" [more ...]" "set timefmt \"<format string>\"
show timefmt") ("timestamp" "set timestamp {\"<format>\"} {top|bottom} {{no}rotate} [more ...]" "set timestamp {\"<format>\"} {top|bottom} {{no}rotate}
              {offset <xoff>{,<yoff>}} {font \"<fontspec>\"}
unset timestamp
show timestamp") ("tics" "set tics {axis | border} {{no}mirror} [more ...]" "set tics {axis | border} {{no}mirror}
         {in | out} {scale {default | <major> {,<minor>}}}
         {{no}rotate {by <ang>}} {offset <offset> | nooffset}
         {left | right | center | autojustify}
         { format \"formatstring\" } { font \"name{,<size>}\" }
         { textcolor <colorspec> }
set tics {front | back}
unset tics
show tics") ("terminal" "set terminal {<terminal-type> | push | pop} [more ...]" "set terminal {<terminal-type> | push | pop}
show terminal") ("table" "set table {\"outfile\"} [more ...]" "set table {\"outfile\"}
plot <whatever>
unset table") ("surface" "set surface [more ...]" "set surface
unset surface
show surface") ("set_style_ellipse" "set style ellipse {units xx|xy|yy} [more ...]" "set style ellipse {units xx|xy|yy}
                  {size {graph|screen} <a>, {{graph|screen} <b>}}
                  {angle <angle>}") ("set_style_rectangle" "set style rectangle {front|back} {lw|linewidth <lw>} [more ...]" "set style rectangle {front|back} {lw|linewidth <lw>}
                    {fillcolor <colorspec>} {fs <fillstyle>}") ("set_style_circle" "set style circle {radius {graph|screen} <R>}") ("set_style_line" "set style line <index> default [more ...]" "set style line <index> default
set style line <index> {{linetype  | lt} <line_type> | <colorspec>}
                       {{linecolor | lc} <colorspec>}
                       {{linewidth | lw} <line_width>}
                       {{pointtype | pt} <point_type>}
                       {{pointsize | ps} <point_size>}
                       {{pointinterval | pi} <interval>}
                       {palette}
unset style line
show style line") ("set_style_increment" "set style increment {default|userstyles} [more ...]" "set style increment {default|userstyles}
show style increment") ("set_style_function" "set style function <plotting-style> [more ...]" "set style function <plotting-style>
show style function") ("set_style_fill" "set style fill {empty [more ...]" "set style fill {empty
                | {transparent} solid {<density>}
                | {transparent} pattern {<n>}}
               {border {lt} {lc <colorspec>} | noborder}") ("set_style_data" "set style data <plotting-style> [more ...]" "set style data <plotting-style>
show style data") ("boxplot_" "set style boxplot {range <r> | fraction <f>} [more ...]" "set style boxplot {range <r> | fraction <f>}
                  {{no}outliers} {pointtype <p>}
                  {candlesticks | financebars}
                  {separation <x>}
                  {labels off | auto | x | x2}
                  {sorted | unsorted}") ("set_style_arrow" "set style arrow <index> default [more ...]" "set style arrow <index> default
set style arrow <index> {nohead | head | heads}
                        {size <length>,<angle>{,<backangle>}}
                        {filled | empty | nofilled}
                        {front | back}
                        { {linestyle | ls <line_style>}
                          | {linetype | lt <line_type>}
                            {linewidth | lw <line_width} }
unset style arrow
show style arrow") ("style" "set style rectangle <object options> <linestyle> <fillstyle> [more ...]" "set style rectangle <object options> <linestyle> <fillstyle>
set style circle radius <size>
set style ellipse size <size> units {xy|xx|yy}") ("style" "set style arrow <n> <arrowstyle> [more ...]" "set style arrow <n> <arrowstyle>
set style fill <fillstyle>
set style histogram <histogram style options>
set style line <n> <linestyle>") ("style" "set style function <style> [more ...]" "set style function <style>
set style data <style>
show style function
show style data") ("size" "set size {{no}square | ratio <r> | noratio} {<xscale>,<yscale>} [more ...]" "set size {{no}square | ratio <r> | noratio} {<xscale>,<yscale>}
show size") ("samples" "set samples <samples_1> {,<samples_2>} [more ...]" "set samples <samples_1> {,<samples_2>}
show samples") ("print_" "set print [more ...]" "set print
set print \"-\"
set print \"<filename>\"
set print \"<filename>\" append
set print \"|<shell_command>\"") ("polar_" "set polar [more ...]" "set polar
unset polar
show polar") ("pointsize" "set pointsize <multiplier> [more ...]" "set pointsize <multiplier>
show pointsize") ("defined_" "set palette  defined { ( <gray1> <color1> {, <grayN> <colorN>}... ) }") ("palette" "set palette [more ...]" "set palette
set palette {
           { gray | color }
           { gamma <gamma> }
           {   rgbformulae <r>,<g>,<b>
             | defined { ( <gray1> <color1> {, <grayN> <colorN>}... ) }
             | file '<filename>' {datafile-modifiers}
             | functions <R>,<G>,<B>
           }
           { cubehelix {start <val>} {cycles <val>} {saturation <val>} }
           { model { RGB | HSV | CMY | YIQ | XYZ } }
           { positive | negative }
           { nops_allcF | ps_allcF }
           { maxcolors <maxcolors> }
         }
show palette
show palette palette <n> {{float | int}}
show palette gradient
show palette fit2rgbformulae
show palette rgbformulae
show colornames") ("parametric_" "set parametric [more ...]" "set parametric
unset parametric
show parametric") ("output" "set output {\"<filename>\"} [more ...]" "set output {\"<filename>\"}
show output") ("origin" "set origin <x-origin>,<y-origin>") ("offsets" "set offsets <left>, <right>, <top>, <bottom> [more ...]" "set offsets <left>, <right>, <top>, <bottom>
unset offsets
show offsets") ("polygon" "set object <index> polygon [more ...]" "set object <index> polygon
    from <position> to <position> ... {to <position>}") ("circle" "set object <index> circle {at|center} <position> size <radius> [more ...]" "set object <index> circle {at|center} <position> size <radius>
    {arc [<begin>:<end>]}
    {<other-object-properties>}") ("ellipse" "set object <index> ellipse {at|center} <position> size <w>,<h> [more ...]" "set object <index> ellipse {at|center} <position> size <w>,<h>
    {angle <orientation>} {units xy|xx|yy}
    {<other-object-properties>}") ("rectangle" "set object <index> rectangle [more ...]" "set object <index> rectangle
    {from <position> {to|rto} <position> |
     center <position> size <w>,<h> |
     at <position> size <w>,<h>}") ("object" "set object <index> [more ...]" "set object <index>
    <object-type> <object-properties>
    {front|back|behind} {fc|fillcolor <colorspec>} {fs <fillstyle>}
    {default} {lw|linewidth <width>}") ("mxtics" "set mxtics {<freq> | default} [more ...]" "set mxtics {<freq> | default}
unset mxtics
show mxtics") ("multiplot" "set multiplot { layout <rows>,<cols> [more ...]" "set multiplot { layout <rows>,<cols>
                {rowsfirst|columnsfirst} {downwards|upwards}
                {title <page title>}
                {scale <xscale>{,<yscale>}} {offset <xoff>{,<yoff>}}
              }
unset multiplot") ("mouse" "set mouse {doubleclick <ms>} {nodoubleclick} \\ [more ...]" "set mouse {doubleclick <ms>} {nodoubleclick} \\
          {{no}zoomcoordinates} \\
          {noruler | ruler {at x,y}} \\
          {polardistance{deg|tan} | nopolardistance} \\
          {format <string>} \\
          {clipboardformat <int>/<string>} \\
          {mouseformat <int>/<string>} \\
          {{no}labels {\"labeloptions\"}} \\
          {{no}zoomjump} {{no}verbose}
unset mouse") ("margin" "set bmargin {{at screen} <margin>} [more ...]" "set bmargin {{at screen} <margin>}
set lmargin {{at screen} <margin>}
set rmargin {{at screen} <margin>}
set tmargin {{at screen} <margin>}
show margin") ("mapping" "set mapping {cartesian | spherical | cylindrical}") ("macros" "set macros") ("logscale" "set logscale <axes> {<base>} [more ...]" "set logscale <axes> {<base>}
unset logscale <axes>
show logscale") ("locale" "set locale {\"<locale>\"}") ("loadpath" "set loadpath {\"pathlist1\" {\"pathlist2\"...}} [more ...]" "set loadpath {\"pathlist1\" {\"pathlist2\"...}}
show loadpath") ("label" "set label {<tag>} {\"<label text>\"} {at <position>} [more ...]" "set label {<tag>} {\"<label text>\"} {at <position>}
          {left | center | right}
          {norotate | rotate {by <degrees>}}
          {font \"<name>{,<size>}\"}
          {noenhanced}
          {front | back}
          {textcolor <colorspec>}
          {point <pointstyle> | nopoint}
          {offset <offset>}
unset label {<tag>}
show label") ("key" "set key {on|off} {default} [more ...]" "set key {on|off} {default}
        {{inside | outside} | {lmargin | rmargin | tmargin | bmargin}
          | {at <position>}}
        {left | right | center} {top | bottom | center}
        {vertical | horizontal} {Left | Right}
        {{no}opaque}
        {{no}reverse} {{no}invert}
        {samplen <sample_length>} {spacing <vertical_spacing>}
        {width <width_increment>}
        {height <height_increment>}
        {{no}autotitle {columnheader}}
        {title \"<text>\"} {{no}enhanced}
        {font \"<face>,<size>\"} {textcolor <colorspec>}
        {{no}box { {linestyle | ls <line_style>}
                   | {linetype | lt <line_type>}
                     {linewidth | lw <line_width>}}}
        {maxcols {<max no. of columns> | auto}}
        {maxrows {<max no. of rows> | auto}}
unset key
show key") ("isosamples" "set isosamples <iso_1> {,<iso_2>} [more ...]" "set isosamples <iso_1> {,<iso_2>}
show isosamples") ("historysize" "set historysize <int> [more ...]" "set historysize <int>
unset historysize") ("hidden3d" "set hidden3d {defaults} | [more ...]" "set hidden3d {defaults} |
             { {front|back}
               {{offset <offset>} | {nooffset}}
               {trianglepattern <bitpattern>}
               {{undefined <level>} | {noundefined}}
               {{no}altdiagonal}
               {{no}bentover} }
unset hidden3d
show hidden3d") ("grid" "set grid {{no}{m}xtics} {{no}{m}ytics} {{no}{m}ztics} [more ...]" "set grid {{no}{m}xtics} {{no}{m}ytics} {{no}{m}ztics}
         {{no}{m}x2tics} {{no}{m}y2tics}
         {{no}{m}cbtics}
         {polar {<angle>}}
         {layerdefault | front | back}
         { {linestyle <major_linestyle>}
           | {linetype | lt <major_linetype>}
             {linewidth | lw <major_linewidth>}
           { , {linestyle | ls <minor_linestyle>}
               | {linetype | lt <minor_linetype>}
                 {linewidth | lw <minor_linewidth>} } }
unset grid
show grid") ("functions_" "show functions") ("format_" "set format {<axes>} {\"<format-string>\"} [more ...]" "set format {<axes>} {\"<format-string>\"}
set format {<axes>} {'<format-string>'}
show format") ("fontpath" "set fontpath {\"pathlist1\" {\"pathlist2\"...}} [more ...]" "set fontpath {\"pathlist1\" {\"pathlist2\"...}}
show fontpath") ("fit_" "set fit {logfile {\"<filename>\"}} {{no}errorvariables} {{no}quiet} [more ...]" "set fit {logfile {\"<filename>\"}} {{no}errorvariables} {{no}quiet}
unset fit
show fit") ("encoding" "set encoding {<value>} [more ...]" "set encoding {<value>}
set encoding locale
show encoding") ("dummy" "set dummy {<dummy-var>} {,<dummy-var>} [more ...]" "set dummy {<dummy-var>} {,<dummy-var>}
show dummy") ("dgrid3d" "set dgrid3d {<rows>} {,{<cols>}} [more ...]" "set dgrid3d {<rows>} {,{<cols>}}
            { splines |
              qnorm {<norm>} |
              (gauss | cauchy | exp | box | hann) 
                {kdensity} {<dx>} {,<dy>} }
unset dgrid3d
show dgrid3d") ("decimalsign" "set decimalsign {<value> | locale {\"<locale>\"}} [more ...]" "set decimalsign {<value> | locale {\"<locale>\"}}
unset decimalsign
show decimalsign") ("set_datafile_binary" "set datafile binary <binary list> [more ...]" "set datafile binary <binary list>
show datafile binary
show datafile
unset datafile") ("set_datafile_commentschars" "set datafile commentschars {\"<string>\"} [more ...]" "set datafile commentschars {\"<string>\"}
show datafile commentschars
unset commentschars") ("set_datafile_separator" "set datafile separator {\"<char>\" | whitespace}") ("set_datafile_missing" "set datafile missing {\"<string>\"} [more ...]" "set datafile missing {\"<string>\"}
show datafile missing
unset datafile") ("contour" "set contour {base | surface | both} [more ...]" "set contour {base | surface | both}
unset contour
show contour") ("cntrparam" "set cntrparam { { linear [more ...]" "set cntrparam { { linear
                | cubicspline
                | bspline
                | points <n>
                | order <n>
                | levels { auto {<n>} | <n>
                           | discrete <z1> {,<z2>{,<z3>...}}
                           | incremental <start>, <incr> {,<end>}
                         }
                }
              }
show contour") ("clip" "set clip <clip-type> [more ...]" "set clip <clip-type>
unset clip <clip-type>
show clip") ("clabel" "set clabel {'<format>'} [more ...]" "set clabel {'<format>'}
unset clabel
show clabel") ("boxwidth" "set boxwidth {<width>} {absolute|relative} [more ...]" "set boxwidth {<width>} {absolute|relative}
show boxwidth") ("border" "set border {<integer>} {front | back} {linewidth | lw <line_width>} [more ...]" "set border {<integer>} {front | back} {linewidth | lw <line_width>}
           {{linestyle | ls <line_style>} | {linetype | lt <line_type>}}
unset border
show border") ("bars" "set bars {small | large | fullwidth | <size>} {front | back} [more ...]" "set bars {small | large | fullwidth | <size>} {front | back}
unset bars
show bars") ("autoscale" "set autoscale {<axes>{|min|max|fixmin|fixmax|fix} | fix | keepfix} [more ...]" "set autoscale {<axes>{|min|max|fixmin|fixmax|fix} | fix | keepfix}
unset autoscale {<axes>}
show autoscale") ("arrow" "set arrow {<tag>} {from <position>} {to|rto <position>} [more ...]" "set arrow {<tag>} {from <position>} {to|rto <position>}
          { {arrowstyle | as <arrow_style>}
            | { {nohead | head | backhead | heads}
                {size <length>,<angle>{,<backangle>}}
                {filled | empty | nofilled}
                {front | back}
                { {linestyle | ls <line_style>}
                  | {linetype | lt <line_type>}
                    {linewidth | lw <line_width} } } }") ("angles" "set angles {degrees | radians} [more ...]" "set angles {degrees | radians}
show angles") ("save" "save  {<option>} '<filename>'") ("raise" "raise {plot_window_nb}") ("print" "print <expression> {, <expression>, ...}") ("with" "with <style> { {linestyle | ls <line_style>} [more ...]" "with <style> { {linestyle | ls <line_style>}
               | {{linetype  | lt <line_type>}
                  {linewidth | lw <line_width>}
                  {linecolor | lc <colorspec>}
                  {pointtype | pt <point_type>}
                  {pointsize | ps <point_size>}
                  {fill | fs <fillstyle>}
                  {nohidden3d} {nocontours} {nosurface}
                  {palette}}
             }") ("title" "title <text> | notitle [<ignored text>] [more ...]" "title <text> | notitle [<ignored text>]
title columnheader | title columnheader(N)") ("iteration_" "plot for [<variable> = <start> : <end> {:<increment>}] [more ...]" "plot for [<variable> = <start> : <end> {:<increment>}]
plot for [<variable> in \"string of words\"]") ("ranges" "[{<dummy-var>=}{{<min>}:{<max>}}] [more ...]" "[{<dummy-var>=}{{<min>}:{<max>}}]
[{{<min>}:{<max>}}]") ("using" "plot 'file' using <entry> {:<entry> {:<entry> ...}} {'format'}") ("thru" "plot 'file' thru f(x)") ("smooth" "smooth {unique | frequency | cumulative | cnormal | kdensity  [more ...]" "smooth {unique | frequency | cumulative | cnormal | kdensity 
               | csplines | acsplines | bezier | sbezier}") ("index" "plot 'file' index { <m>{:<n>{:<p>}} | \"<name>\" }") ("every" "plot 'file' every {<point_incr>} [more ...]" "plot 'file' every {<point_incr>}
                    {:{<block_incr>}
                      {:{<start_point>}
                        {:{<start_block>}
                          {:{<end_point>}
                            {:<end_block>}}}}}") ("data" "plot '<file_name>' {binary <binary list>} [more ...]" "plot '<file_name>' {binary <binary list>}
                   {{nonuniform} matrix}
                   {index <index list> | index \"<name>\"}
                   {every <every list>}
                   {thru <thru expression>}
                   {using <using list>}
                   {smooth <option>}
                   {volatile} {noautoscale}") ("general" "plot '<file_name>' {binary <binary list>} ... [more ...]" "plot '<file_name>' {binary <binary list>} ...
splot '<file_name>' {binary <binary list>} ...") ("plot" "plot {<ranges>} [more ...]" "plot {<ranges>}
     {<iteration>}
     {<function> | {\"<datafile>\" {datafile-modifiers}}}
     {axes <axes>} {<title-spec>} {with <style>}
     {, {definitions{,}} <function> ...}") ("pause" "pause <time> {\"<string>\"} [more ...]" "pause <time> {\"<string>\"}
pause mouse {<endcondition>}{, <endcondition>} {\"<string>\"}") ("lower" "lower {plot_window_nb}") ("load" "load \"<input-file>\"") ("fit" "fit {<ranges>} <expression> [more ...]" "fit {<ranges>} <expression>
    '<datafile>' {datafile-modifiers}
    via '<parameter file>' | <var1>{,<var2>,...}") ("evaluate" "eval <string expression>") ("Do" "do for <iteration-spec> { [more ...]" "do for <iteration-spec> {
     <commands>
     <commands>
}") ("call" "call \"<input-file>\" <parameter-0> <parm-1> ... <parm-9>") ("cd" "cd '<directory-name>'") ("newhistogram" "newhistogram {\"<title>\"} {lt <linetype>} {fs <fillstyle>} {at <x-coord>}") ("filledcurves" "plot ... with filledcurves [option]") ("bind" "bind {allwindows} [<key-sequence>] [\"<gnuplot commands>\"] [more ...]" "bind {allwindows} [<key-sequence>] [\"<gnuplot commands>\"]
bind <key-sequence> \"\"
reset bind") ("colorspec" "... {linecolor | lc} {<colorspec> | <n>} [more ...]" "... {linecolor | lc} {<colorspec> | <n>}
... {textcolor | tc} {<colorspec> | {linetype | lt} <n>}"))))) (while alist (puthash (caar alist) (cdar alist) tbl) (setq alist (cdr alist))) tbl))