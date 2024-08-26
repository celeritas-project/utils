#include <array>
#include <cmath>
#include <iostream>
#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTree.h>

#include "../src/RootData.hh"

//---------------------------------------------------------------------------//
/*!
 * cos(theta) of two normalized vectors: cos(a,b) = dot(a,b) / |a||b|
 */
double cos_theta(rootdata::Array3 const& a, rootdata::Array3 const& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;  // |a| = |b| = 1
}
//---------------------------------------------------------------------------//
void verify_pol()
{
    gSystem->Load("librootdata");

    auto file = TFile::Open("out-lar.1e-5mev.root", "read");
    auto tree = file->Get<TTree>("events");

    rootdata::Event* event = nullptr;
    tree->SetBranchAddress("event", &event);

    auto h = new TH1D("", "", 100, 0, TMath::Pi());

    for (auto i = 0; i < tree->GetEntries(); i++)
    {
        tree->GetEntry(i);
        auto const& p = event->primaries[0];
        auto const& inc = p.steps[0].polarization;
        auto const& out = p.steps[1].polarization;

        auto result = cos_theta(inc, out);
        if (result < 0.999999)
        {
            h->Fill(std::acos(result));
        }
    }

    auto c = new TCanvas("", "", 700, 500);
    h->Draw();
}
