#ifndef __KMAPWidgetWorkers_h
#define __KMAPWidgetWorkers_h
#include <QThread>
#include "vtkSlicerKMAPLogic.h"
#include <QString>


class TCMWorker : public QThread
{
    Q_OBJECT
public:
    TCMWorker(vtkSlicerKMAPLogic* logic,
              std::vector<std::vector<double>>& voxels,
              std::vector<double>& Cp,
              std::vector<double>& framing,
              const std::vector<std::string>& modelIDs,
              const double& vbInit, const double& vbLower, double& vbUpper,
              const double& k1Init, const double& k1Lower, double& k1Upper,
              const double& k2Init, const double& k2Lower, double& k2Upper,
              const double& k3Init, const double& k3Lower, double& k3Upper,
              const double& k4Init, const double& k4Lower, double& k4Upper,
              const double& tdInit, const double& tdLower, double& tdUpper,
              double dk,
              double timestep,
              const double pbrp[3],
              int maxiter,
              std::atomic<bool>& stopRequested,
              const std::vector<double>* wgt_global = nullptr,
              int numThreads = 1)
        : logic(logic),
          voxels(voxels),
          Cp(Cp),
          framing(framing),
          modelIDs(modelIDs),
          dk(dk),
          timestep(timestep),
          maxiter(maxiter),
          stopRequested(stopRequested),
          numThreads(numThreads)
    {
      std::copy(pbrp, pbrp+3, this->pbrp);
      if (wgt_global) wgt_copy = *wgt_global;
      // copy init/lb/ub arrays
      this->vbInit = vbInit; this->vbLower = vbLower; this->vbUpper = vbUpper;
      this->k1Init = k1Init; this->k1Lower = k1Lower; this->k1Upper = k1Upper;
      this->k2Init = k2Init; this->k2Lower = k2Lower; this->k2Upper = k2Upper;
      this->k3Init = k3Init; this->k3Lower = k3Lower; this->k3Upper = k3Upper;
      this->k4Init = k4Init; this->k4Lower = k4Lower; this->k4Upper = k4Upper;
      this->tdInit = tdInit; this->tdLower = tdLower; this->tdUpper = tdUpper;
    }

    std::vector<TCMParameters> results;

signals:
    void progressChanged(int value);
    void finishedProcessing(const QString& modelID);//, const std::vector<TCMParameters>& results);
    void canceled(const QString& modelID);
    void finishedAll();
    void modelStarted(const QString& modelID);

protected:
    void run() override
    {
        for (size_t i = 0; i < modelIDs.size(); ++i) {
          if (stopRequested) break;

          const std::string& modelID = modelIDs[i];
          emit modelStarted(QString::fromStdString(modelID));

          // 1. Configure model parameters
          bool sens[6] = {false};
          double init[6]{0}, lb[6]{0}, ub[6]{0};
          int num_tc = 0;

          configureModel(modelID, sens, init, lb, ub, num_tc, modelfields);

          std::vector<TCMParameters> localOutput;
          logic->callTCMImg(voxels,
                            Cp,
                            framing,
                            init,
                            lb,
                            ub,
                            sens,
                            dk,
                            timestep,
                            pbrp,
                            maxiter,
                            num_tc,
                            localOutput,
                            modelID,
                            stopRequested,
                            &wgt_copy,
                            numThreads,
                            [this](int value){ emit progressChanged(value); },
                            [this](){ return stopRequested.load(); });
            if (stopRequested) {
                emit canceled(QString::fromStdString(modelID));
                break;
            } else {
                results = std::move(localOutput);
                emit finishedProcessing(QString::fromStdString(modelID));
            }
        }
        emit finishedAll();
    }

private:
    void configureModel(const std::string& modelID, bool sens[6], double init[6], double lb[6], double ub[6], int& num_tc, std::vector<std::string>& fields)
    {
      if (modelID == "1TCM") {
          num_tc = 1;
          sens[0] = true; sens[1] = true; sens[2] = true; sens[3] = false;
          lb[0] = vbLower; lb[1] = k1Lower; lb[2] = k2Lower;
          ub[0] = vbUpper; ub[1] = k1Upper; ub[2] = k2Upper;
          init[0] = vbInit; init[1] = k1Init; init[2] = k2Init;
          modelfields = {"K1", "k2", "vb", "Ki", "DV", "AIC", "MASE", "BIC", "chi2"};
      }
      else if (modelID == "1TdCM") {
          num_tc = 1;
          sens[0] = true; sens[1] = true; sens[2] = true; sens[3] = true;
          lb[0] = vbLower; lb[1] = k1Lower; lb[2] = k2Lower; lb[3] = tdLower;
          ub[0] = vbUpper; ub[1] = k1Upper; ub[2] = k2Upper; ub[3] = tdUpper;
          init[0] = vbInit; init[1] = k1Init; init[2] = k2Init; init[3] = tdInit;
          modelfields = {"K1", "k2", "vb", "td", "Ki", "DV", "AIC", "MASE", "BIC", "chi2"};
      }
      else if (modelID == "1TiCM") {
          num_tc = 1;
          sens[0] = true; sens[1] = true;
          lb[0] = vbLower; lb[1] = k1Lower;
          ub[0] = vbUpper; ub[1] = k1Upper;
          init[0] = vbInit; init[1] = k1Init;
          modelfields = {"K1", "vb", "Ki", "AIC", "MASE", "BIC", "chi2"};
      }
      else if (modelID == "1TidCM") {
          num_tc = 1;
          sens[0] = true; sens[1] = true; sens[3] = true;
          lb[0] = vbLower; lb[1] = k1Lower; lb[3] = tdLower;
          ub[0] = vbUpper; ub[1] = k1Upper; ub[3] = tdUpper;
          init[0] = vbInit; init[1] = k1Init; init[3] = tdInit;
          modelfields = {"K1", "vb", "td", "Ki", "AIC", "MASE", "BIC", "chi2"};
      }
      else if (modelID == "2TCM") {
          num_tc = 2;
          sens[0] = true; sens[1] = true; sens[2] = true; sens[3] = true; sens[4] = true;
          lb[0] = vbLower; lb[1] = k1Lower; lb[2] = k2Lower; lb[3] = k3Lower; lb[4] = k4Lower;
          ub[0] = vbUpper; ub[1] = k1Upper; ub[2] = k2Upper; ub[3] = k3Upper; ub[4] = k4Upper;
          init[0] = vbInit; init[1] = k1Init; init[2] = k2Init; init[3] = k3Init; init[4] = k4Init;
          modelfields = {"K1", "k2", "k3", "k4", "vb", "Ki", "DV", "AIC", "MASE", "BIC", "chi2"};
      }
      else if (modelID == "2dTCM") {
          num_tc = 2;
          sens[0] = true; sens[1] = true; sens[2] = true; sens[3] = true; sens[4] = true; sens[5] = true;
          lb[0] = vbLower; lb[1] = k1Lower; lb[2] = k2Lower; lb[3] = k3Lower; lb[4] = k4Lower; lb[5] = tdLower;
          ub[0] = vbUpper; ub[1] = k1Upper; ub[2] = k2Upper; ub[3] = k3Upper; ub[4] = k4Upper; ub[5] = tdUpper;
          init[0] = vbInit; init[1] = k1Init; init[2] = k2Init; init[3] = k3Init; init[4] = k4Init; init[5] = tdInit;
          modelfields = {"K1", "k2", "k3", "k4", "vb", "td", "Ki", "DV", "AIC", "MASE", "BIC", "chi2"};
      }
      else if (modelID == "2TiCM") {
          num_tc = 2;
          sens[0] = true; sens[1] = true; sens[2] = true; sens[3] = true;
          lb[0] = vbLower; lb[1] = k1Lower; lb[2] = k2Lower; lb[3] = k3Lower;
          ub[0] = vbUpper; ub[1] = k1Upper; ub[2] = k2Upper; ub[3] = k3Upper;
          init[0] = vbInit; init[1] = k1Init; init[2] = k2Init; init[3] = k3Init;
          modelfields = {"K1", "k2", "k3", "vb", "Ki", "AIC", "MASE", "BIC", "chi2"};
      }
      else if (modelID == "2TidCM") {
          num_tc = 2;
          sens[0] = true; sens[1] = true; sens[2] = true; sens[3] = true; sens[5] = true;
          lb[0] = vbLower; lb[1] = k1Lower; lb[2] = k2Lower; lb[3] = k3Lower; lb[5] = tdLower;
          ub[0] = vbUpper; ub[1] = k1Upper; ub[2] = k2Upper; ub[3] = k3Upper; ub[5] = tdUpper;
          init[0] = vbInit; init[1] = k1Init; init[2] = k2Init; init[3] = k3Init; init[5] = tdInit;
          modelfields = {"K1", "k2", "k3", "vb", "td", "Ki", "AIC", "MASE", "BIC", "chi2"};
      }
      else {
          throw std::runtime_error(
              "Unknown model ID " + modelID
          );
          return;
      }
    }

    vtkSlicerKMAPLogic* logic;
    std::vector<std::vector<double>> voxels;
    std::vector<double> Cp;
    std::vector<double> framing;
    std::vector<std::string> modelIDs;

    std::atomic<bool>& stopRequested;
    double dk, timestep;
    double pbrp[3];
    int maxiter, numThreads;
    std::vector<double> wgt_copy;

    // All init/lb/ub vectors
    std::vector<std::string> modelfields;
    double vbInit, vbLower, vbUpper;
    double k1Init, k1Lower, k1Upper;
    double k2Init, k2Lower, k2Upper;
    double k3Init, k3Lower, k3Upper;
    double k4Init, k4Lower, k4Upper;
    double tdInit, tdLower, tdUpper;
};

#endif
