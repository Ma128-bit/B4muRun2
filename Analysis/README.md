# Analysis code for the B&rarr;4mu search at Run 2

## Run analysis on an year with condor (Data and MC):
```
source prepare_condor.sh [year] [Analysis_type] [delta]
```
*  [year] is the year (`2017` or `2018`);
*  [Analysis_type] `B4mu`, `Norm`;
*  [Delta] is the number of input files per submission

**FOR DATA ONLY**:
```
cd [Analysis_type]/[year]
source submit_era.sh
```
Then:
```
source hadd_era.sh
```
**FOR MC**:
```
cd [Analysis_type]/[year]_[era]/[MC_lable]
submit submit.condor with condor
merge the output with hadd
```
<p>&nbsp;</p>

## Run analysis on a specific subset of data

```
python3 Analizer.py --index [ID] --delta [Delta] --directory_IN [Input_dir] --directory_OUT [Output_dir] --isMC [isMC] --Analysis_type [Analysis_type]
```
*  [Input_dir] is the directory with the root files;
*  [Output_dir] is the directory where the files will be saved;
*  [Delta] is the number of input files;
*  [ID] is used to select input file root with number between [ID]x[Delta] and ([ID]+1)x[Delta];
*  [isMC] is 0 for data and 1 for MC;
*  [Analysis_type] `B4mu`, `Norm`

<p>&nbsp;</p>
