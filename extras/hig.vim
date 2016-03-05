" VIM syntax file
" Language: HipGISAXS input file format
" Maintainer: Abhinav Sarje <asarje@lbl.gov>
" Created: June 09 2012
" Revised: March 05 2016

if exists("b:current_syntax")
	finish
endif

" matches
syn match higNumber '\d\+'
syn match higNumber '[-+]\d\+'
syn match higNumber '\d\+\.\d*'
syn match higNumber '[-+]\d\+\.\d*'
syn match higNumber '[-+]\=\d[[:digit:]]*[eE][\-+]\=\d\+'
syn match higNumber '\d[[:digit:]]*[eE][\-+]\=\d\+'
syn match higNumber '[-+]\=\d[[:digit:]]*\.\d*[eE][\-+]\=\d\+'
syn match higNumber '\d[[:digit:]]*\.\d*[eE][\-+]\=\d\+'

"syntax region higQuoted matchgroup=higQuotedDelimiter start=/[^\[\]|\\ ]\+\[=/ end=/=\]/ skip=/\\\[=\|\\\=\]/ contains=higEscape,higParamSeparator
syntax match higParamSeparator /,/
syntax match higEscape /\\./

syn region    higString     start=+L\="+ skip=+\\\\\|\\"+ end=+"+ contains=@Spell

" For comments
syn region    higComment   start="#" skip="\\$" end="$" keepend contains=@Spell

" language keywords
syn keyword higMain hipGisaxsInput
syn keyword higMainComponents shape layer unitcell structure instrumentation computation fitting
syn keyword higShapeComponents key name originvec ztilt xyrotation param
syn keyword higLayerComponents key order thickness refindex
syn keyword higUnitcellComponents key element locations
syn keyword higStructureComponents key grain ensemble
syn keyword higGrainComponents shape:key layer:key refindex lattice scaling transvec repetition dimensions xspacing yspacing domain volfraction
syn keyword higEnsembleComponents spacing maxgrains distribution orientations
syn keyword higInstrumentationComponents scattering detector
syn keyword higScatteringComponents expt alphai inplanerot tilt photon polarization coherence spotarea smearing
syn keyword higDetectorComponents origin totalpixels pixelsize sdd directbeam
syn keyword higComputationComponents pathprefix inputdir runname method outputregion resolution nslices structcorrelation saveff savesf
syn keyword higFittingComponenets fitparam key variable range init referencedata algorithm path fitregion npoints algoname algoorder algoparam restart tolerance
syn keyword higShapeParam type min max stat p1 p2 nvalues nextgroup=higNumber skipwhite
syn keyword higRefindexParam delta beta
syn keyword higLatticeParam a b c type hkl abangle caratio
syn keyword higRotationParam axis angles
syn keyword higOrientationParam rot1 rot2 rot3
syn keyword higTiltParam min max step
syn keyword higPhotonParam value unit
syn keyword higOutputParam type minpoint maxpoint
syn keyword higDistributionParam mean stddev
syn keyword higFittingParam type regmin regmax parallel perpendicular pvalue

let b:current_syntax = "hig"
hi def link higNumber						Number
hi def link higMain							Structure
hi def link higMainComponents				Function
hi def link higShapeComponents				Function
hi def link higLayerComponents				Function
hi def link higUnitcellComponents				Function
hi def link higStructureComponents			Function
hi def link higGrainComponents				Function
hi def link higEnsembleComponents			Function
hi def link higInstrumentationComponents	Function
hi def link higScatteringComponents			Function
hi def link higDetectorComponents			Function
hi def link higComputationComponents		Function
hi def link higFittingComponenets			Function
hi def link higShapeParam					Type
hi def link higRefindexParam				Type
hi def link higLatticeParam					Type
hi def link higOrientationParam				Type
hi def link higRotationParam				Type
hi def link higTiltParam					Type
hi def link higPhotonParam					Type
hi def link higOutputParam					Type
hi def link higDistributionParam    Type
hi def link higFittingParam					Type
"highlight link higQuoted					String
hi def link higEscape						Special
hi def link higParamSeparator				Delimiter
hi def link higQuotedDelimiter				Special
hi def link higString						String
hi def link higComment						Comment
