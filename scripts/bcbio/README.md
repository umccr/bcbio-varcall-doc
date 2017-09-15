* `2017-08-28_cancer-normal`: bcbio run with VarDict
    * Ran smoothly
* `2017-09-11_cancer-normal`: bcbio run with MuTect2.
    * Needed to add `tools_on: [gatk4]` in the config file.
    * Needed to add `tools_off: [gemini]` in the config file.
