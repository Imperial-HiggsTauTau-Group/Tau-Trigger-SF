/*! Definition of c++ methods used in python code.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#pragma once

enum class LegType { e = 1, mu = 2, tau = 4, jet = 8 };

class PileUpWeightProvider {
public:
    PileUpWeightProvider(const TH1D& data_pu, const TH1D& mc_pu)
    {
        ratio.reset(static_cast<TH1D*>(mc_pu.Clone("pu_weight_ratio")));
        ratio->SetDirectory(nullptr);
        ratio->Reset();
        const double data_norm = data_pu.Integral(), mc_norm = mc_pu.Integral();
        for (int i = 1; i <= mc_pu.GetNbinsX(); ++i) {
            const double lo = mc_pu.GetXaxis()->GetBinLowEdge(i), hi = mc_pu.GetXaxis()->GetBinUpEdge(i);
            const int j = data_pu.GetXaxis()->FindFixBin(0.5 * (lo + hi));
            if (std::abs(data_pu.GetXaxis()->GetBinLowEdge(j) - lo) > 1e-6
                    || std::abs(data_pu.GetXaxis()->GetBinUpEdge(j) - hi) > 1e-6)
                throw std::runtime_error("PileUpWeightProvider: data and MC pileup bin edges do not align.");
            const double p_mc = mc_pu.GetBinContent(i) / mc_norm;
            ratio->SetBinContent(i, p_mc > 0 ? data_pu.GetBinContent(j) / data_norm / p_mc : 0.);
        }
    }

    float GetWeight(float npu) const
    {
        const int bin = ratio->FindFixBin(npu);
        if (bin < 1 || bin > ratio->GetNbinsX())
            return 0;
        return ratio->GetBinContent(bin);
    }

    static void Initialize(const TH1D& data_pu, const TH1D& mc_pu)
    {
        default_provider.reset(new PileUpWeightProvider(data_pu, mc_pu));
    }

    static const PileUpWeightProvider& GetDefault()
    {
        if(!default_provider)
            throw std::runtime_error("Default PileUpWeightProvider is not initialized.");
        return *default_provider;
    }

private:
    static std::unique_ptr<PileUpWeightProvider> default_provider;

private:
    std::unique_ptr<TH1D> ratio;
};

std::unique_ptr<PileUpWeightProvider> PileUpWeightProvider::default_provider;

