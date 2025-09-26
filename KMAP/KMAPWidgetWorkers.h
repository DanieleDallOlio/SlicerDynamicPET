#ifndef __KMAPWidgetWorkers_h
#define __KMAPWidgetWorkers_h
#include <QThread>
#include "vtkSlicerKMAPLogic.h"


class TCMWorker : public QThread
{
    Q_OBJECT
public:
    TCMWorker(vtkSlicerKMAPLogic* logic,
              std::vector<std::vector<double>> voxels,
              std::vector<double> Cp,
              std::vector<double> framing,
              const double* kinit,
              const double* lb,
              const double* ub,
              const bool* sens,
              double dk,
              double timestep,
              const double pbrp[],
              int maxiter,
              int n_tc,
              const QString& modelID,
              std::atomic<bool>& stopRequested,
              const std::vector<double>* wgt_global,
              int numThreads)
        : logic(logic),
          voxels(std::move(voxels)),
          Cp(std::move(Cp)),
          framing(std::move(framing)),
          dk(dk),
          timestep(timestep),
          maxiter(maxiter),
          n_tc(n_tc),
          modelID(modelID),
          stopRequested(stopRequested),
          numThreads(numThreads)
    {
        // number of parameters depends on model
        int num_par = (n_tc == 1) ? 4 : 6;

        // copy arrays so they stay valid
        kinit_copy.assign(kinit, kinit + num_par);
        lb_copy.assign(lb, lb + num_par);
        ub_copy.assign(ub, ub + num_par);

        sens_copy.resize(num_par);
        for (int i = 0; i < num_par; ++i)
            sens_copy[i] = sens[i] ? 1 : 0;

        // pbrp has length 3 in your code
        pbrp_copy.assign(pbrp, pbrp + 3);

        if (wgt_global) {
            wgt_copy = *wgt_global;  // deep copy of weights
            wgt_copy_ptr = &wgt_copy; // safe pointer for call
        } else {
            wgt_copy_ptr = nullptr;
        }
    }

protected:
    void run() override
    {
        std::vector<TCMParameters> localOutput;
        logic->callTCMImg(voxels,
                          Cp,
                          framing,
                          kinit_copy.data(),
                          lb_copy.data(),
                          ub_copy.data(),
                          reinterpret_cast<const bool*>(sens_copy.data()),
                          dk,
                          timestep,
                          pbrp_copy.data(),
                          maxiter,
                          n_tc,
                          localOutput,
                          modelID.toStdString(),
                          stopRequested,
                          wgt_copy_ptr,
                          numThreads,
                          [this](int value){ emit progressChanged(value); },
                          [this](){ return stopRequested.load(); });

        if (stopRequested) {
          emit canceled(modelID);
        } else {
          emit finishedProcessing(modelID, localOutput);
        }
    }

signals:
    void progressChanged(int value);
    void finishedProcessing(const QString& modelID,
                            const std::vector<TCMParameters>& results);
    void canceled(const QString& modelID);

private:
    vtkSlicerKMAPLogic* logic;
    std::vector<std::vector<double>> voxels;
    std::vector<double> Cp;
    std::vector<double> framing;
    double dk;
    double timestep;
    int maxiter;
    int n_tc;
    QString modelID;
    std::atomic<bool>& stopRequested;
    const std::vector<double>* wgt_global_copy;
    int numThreads;

    // safe copies of arrays
    std::vector<double> kinit_copy, lb_copy, ub_copy, pbrp_copy;
    std::vector<char> sens_copy; // avoids vector<bool> weirdness
    std::vector<double> wgt_copy;
    const std::vector<double>* wgt_copy_ptr = nullptr;
};

#endif
