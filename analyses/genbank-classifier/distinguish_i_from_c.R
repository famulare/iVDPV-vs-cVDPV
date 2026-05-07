# Build the reference iVDPV-vs-cVDPV classifier and regenerate its reusable
# plots/tables/artifacts. The implementation lives in utils/R/classifier_workflow.R so
# follow-up explorations can reuse individual steps without copying this script.
source(file.path("utils", "R", "classifier_workflow.R"))

run_classifier_workflow(classifier_config())
