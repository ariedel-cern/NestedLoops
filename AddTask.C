/**
 * File              : AddTask.C
 * Author            : Anton Riedel <anton.riedel@tum.de>
 * Date              : 07.05.2021
 * Last Modified Date: 23.11.2021
 * Last Modified By  : Anton Riedel <anton.riedel@tum.de>
 */

// AddTask* macro for flow analysis
// Creates a Flow Event task and adds it to the analysis manager.
// Sets the cuts using the correction framework (CORRFW) classes.
// Also creates Flow Analysis tasks and connects them to the output of the flow
// event task.

void AddTask(Float_t centerMin = 0., Float_t centerMax = 100.,
             Bool_t bRunOverAOD = kTRUE) {

  TString OutputFile(std::getenv("GRID_OUTPUT_ROOT_FILE"));

  // Get the pointer to the existing analysis manager
  // via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask.C macro", "No analysis manager to connect to.");
    return;
  }

  // Check the analysis type using the event handlers connected to the analysis
  // manager. The availability of MC handler can also be checked here.
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask.C macro", "This task requires an input event handler");
    return;
  }

  AliAnalysisTaskAR *task =
      new AliAnalysisTaskAR(Form("%s%s_%.1f-%.1f", std::getenv("TASK_BASENAME"),
                                 "Qvector", centerMin, centerMax));

  task->SetDefaultConfiguration();
  task->SetDefaultBinning();
  task->SetDefaultCuts(128, centerMin, centerMax);

  std::vector<Double_t> ptBinning = {0.2, 0.34, 0.5, 0.7, 1.0, 2.0, 5.0};
  task->SetTrackBinning(kPT, ptBinning);
  std::vector<Double_t> etaBinning = {-0.8, -0.4, 0., 0.4, 0.8};
  task->SetTrackBinning(kETA, etaBinning);

  // setters for symmetric cumulants we want to compute
  std::vector<std::vector<Int_t>> correlator = {{-2, 2}};
  task->SetCorrelators(correlator);
  task->SetFixedMultiplicity(static_cast<Int_t>(30));

  // task with nested loops
  AliAnalysisTaskAR *taskNestedLoops = dynamic_cast<AliAnalysisTaskAR *>(
      task->Clone(Form("%s%s_%.1f-%.1f", std::getenv("TASK_BASENAME"),
                       "NestedLoops", centerMin, centerMax)));
  taskNestedLoops->SetUseNestedLoops(kTRUE);

  // generate some weight histograms
  Int_t phi_bins = 100;
  Double_t phi_min = 0;
  Double_t phi_max = TMath::TwoPi();
  Double_t phi2 = TMath::TwoPi() / 6;
  Double_t phi3 = TMath::TwoPi() / 3;
  Double_t phiWeight = 1.4;
  TH1D *WeightPhi =
      new TH1D("phi_weight", "phi_weight", phi_bins, phi_min, phi_max);
  for (int i = 0; i < WeightPhi->GetNbinsX(); i++) {
    if (phi2 < WeightPhi->GetBinCenter(i + 1) &&
        phi3 > WeightPhi->GetBinCenter(i + 1)) {
      WeightPhi->SetBinContent(i + 1, phiWeight);
    } else {
      WeightPhi->SetBinContent(i + 1, 1.);
    }
  }

  Int_t pt_bins = 100;
  Double_t pt_min = 0.2;
  Double_t pt_max = 5.;
  Double_t pt2 = 0.4;
  Double_t pt3 = 1.2;
  Double_t ptWeight = 1.6;
  TH1D *WeightPt = new TH1D("pt_weight", "pt_weight", pt_bins, pt_min, pt_max);
  for (int i = 0; i < WeightPt->GetNbinsX(); i++) {
    if (pt2 < WeightPt->GetBinCenter(i + 1) &&
        pt3 > WeightPt->GetBinCenter(i + 1)) {
      WeightPt->SetBinContent(i + 1, ptWeight);
    } else {
      WeightPt->SetBinContent(i + 1, 1.);
    }
  }

  Int_t eta_bins = 100;
  Double_t eta_min = -0.8;
  Double_t eta_max = 0.8;
  Double_t eta2 = -0.1;
  Double_t eta3 = 0.4;
  Double_t etaWeight = 2.4;
  TH1D *WeightEta =
      new TH1D("eta_weight", "eta_weight", eta_bins, eta_min, eta_max);
  for (int i = 0; i < WeightEta->GetNbinsX(); i++) {
    if (eta2 < WeightEta->GetBinCenter(i + 1) &&
        eta3 > WeightEta->GetBinCenter(i + 1)) {
      WeightEta->SetBinContent(i + 1, etaWeight);
    } else {
      WeightEta->SetBinContent(i + 1, 1.);
    }
  }

  // task with weights
  AliAnalysisTaskAR *taskWithWeights = dynamic_cast<AliAnalysisTaskAR *>(
      task->Clone(Form("%s%s_%.1f-%.1f", std::getenv("TASK_BASENAME"),
                       "QVectorWithWeights", centerMin, centerMax)));
  taskWithWeights->SetWeightHistogram(kPHI, WeightPhi);
  taskWithWeights->SetWeightHistogram(kPT, WeightPt);
  taskWithWeights->SetWeightHistogram(kETA, WeightEta);

  AliAnalysisTaskAR *taskNestedLoopsWithWeights =
      dynamic_cast<AliAnalysisTaskAR *>(taskWithWeights->Clone(
          Form("%s%s_%.1f-%.1f", std::getenv("TASK_BASENAME"),
               "NestedLoopsWithWeights", centerMin, centerMax)));
  taskNestedLoopsWithWeights->SetUseNestedLoops(kTRUE);

  // add all tasks to the analysis manager in a loop
  std::vector<AliAnalysisTaskAR *> tasks = {
      task, taskNestedLoops, taskWithWeights, taskNestedLoopsWithWeights};

  // Define input/output containers:
  OutputFile += TString(":") + TString(std::getenv("OUTPUT_TDIRECTORY_FILE"));
  AliAnalysisDataContainer *cinput = nullptr;
  AliAnalysisDataContainer *coutput = nullptr;

  // loop over all tasks
  for (auto T : tasks) {
    mgr->AddTask(T);
    cout << "Added to manager: " << T->GetName() << endl;
    cinput = mgr->GetCommonInputContainer();
    coutput =
        mgr->CreateContainer(T->GetName(), TList::Class(),
                             AliAnalysisManager::kOutputContainer, OutputFile);
    mgr->ConnectInput(T, 0, cinput);
    mgr->ConnectOutput(T, 1, coutput);
  }

  std::cout << "YEAH" << std::endl;
}
