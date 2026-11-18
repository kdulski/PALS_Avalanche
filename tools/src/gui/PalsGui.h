#ifndef PALS_GUI_H
#define PALS_GUI_H

// Graphical interface of PALS_avalanche (ROOT GUI classes).
//
// Left: the spectrum (log scale) with the fit, components, fit range and background
// window, residuals below. Mouse modes: zoom (ROOT default), select the fit range,
// select the background window. Right: files, data settings, range and estimates,
// resolution and lifetime components, fit settings, fit / save buttons, report.
//
// Events are handled with ProcessMessage (widgets) and a TCanvas subclass (mouse),
// so no ROOT dictionary is needed.

#include "Config.h"
#include "FitSession.h"

#include <TCanvas.h>
#include <TGButton.h>
#include <TGFrame.h>

#include <memory>
#include <string>
#include <vector>

class TBox;
class TGCheckButton;
class TGComboBox;
class TGLabel;
class TGNumberEntry;
class TGRadioButton;
class TGGC;
class TGTextEntry;
class TGTextView;
class TH1D;
class TPad;
class TRootEmbeddedCanvas;

namespace pals {

class PalsGui;

// Text button drawn as a flat rounded rectangle (classic ROOT buttons are square, 3D).
class RoundButton : public TGTextButton {
public:
    RoundButton(const TGWindow* parent, const char* text, Int_t id = -1);
    ~RoundButton() override;
    void ChangeBackground(Pixel_t back) override;   // colour around the rounded corners
    void DrawBorder() override {}
    Bool_t HandleCrossing(Event_t* event) override;
    TGDimension GetDefaultSize() const override;

protected:
    void DoRedraw() override;

private:
    TGGC* gc_ = nullptr;
    Pixel_t surround_ = 0;
    bool hover_ = false;
};

// Colours of the interface (one place to change the look).
struct GuiTheme {
    static Pixel_t panel();          // background of the control area (#87ceeb, sky blue)
    static Pixel_t text();           // labels, buttons, group titles (navy grey)
    static Pixel_t button();         // button face
    static Pixel_t buttonHover();
    static Pixel_t buttonDown();
    static Pixel_t buttonDisabled();
    static Pixel_t buttonBorder();
    static FontStruct_t font(bool bold = false);   // sans-serif (Helvetica or a fallback), 0 = keep ROOT's
    static constexpr int kCornerRadius = 7;
    // where the PAv logo is shown (trial: all three; switch off the ones you do not want)
    static constexpr bool kLogoInFitPanel = false;    // card on the right of the "Fit" settings
    static constexpr bool kLogoBottomLeft = true;     // small mark below the residuals, left of the status
    static constexpr bool kLogoWatermark = true;      // faint logo behind the main plot
    static constexpr int kBorderWidth = 2;
    static constexpr int kButtonExtraWidth = 4;    // pixels added to the default button width
    static constexpr int kButtonExtraHeight = 1;   // and height
};

// Canvas that forwards mouse selections to the GUI.
class SpectrumCanvas : public TCanvas {
public:
    SpectrumCanvas(const char* name, Int_t ww, Int_t wh, Int_t winid, PalsGui* gui);
    void HandleInput(EEventType event, Int_t px, Int_t py) override;

private:
    PalsGui* gui_;
    bool selecting_ = false;
    double x0_ = 0.0;
};

class PalsGui : public TGMainFrame {
public:
    enum MouseMode { kZoom, kFitRange, kBackground };

    PalsGui(const TGWindow* parent, UInt_t w, UInt_t h);
    ~PalsGui() override;

    void openConfig(const std::string& path);
    void openData(const std::string& path);

    // from SpectrumCanvas
    MouseMode mouseMode() const { return mouseMode_; }
    TPad* spectrumPad() const { return padSpectrum_; }
    void rangeSelected(double x1, double x2);
    void selectionMoved(double x1, double x2);
    void cursorMoved(double x);

    Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2) override;
    void CloseWindow() override;

private:
    struct ResolutionRow {
        TGHorizontalFrame* frame;
        TGNumberEntry* sigma;
        TGCheckButton* sigmaAuto;
        TGNumberEntry* intensity;
        TGComboBox* sigmaMode;
        TGComboBox* intensityMode;
    };
    struct LifetimeRow {
        TGHorizontalFrame* frame;
        TGNumberEntry* tau;
        TGNumberEntry* intensity;
        TGComboBox* tauMode;
        TGComboBox* intensityMode;
        TGCheckButton* pPs;
        TGCheckButton* hasWidth;
        TGNumberEntry* width;
        TGComboBox* widthMode;
    };

    void buildLeft(TGCompositeFrame* parent);
    void buildRight(TGCompositeFrame* parent);
    void applyTheme();
    TGGroupFrame* group(TGCompositeFrame* parent, const char* title);

    void populate(const Config& cfg);            // config -> widgets
    Config collect() const;                      // widgets -> config (validated)
    void rebuildResolutionRows(const std::vector<ResolutionComponent>& rows);
    void rebuildLifetimeRows(const std::vector<LifetimeComponent>& rows);
    std::vector<ResolutionComponent> currentResolution() const;
    std::vector<LifetimeComponent> currentLifetimes() const;

    void applyAndEstimate();                     // collect, (re)load data if needed, estimate, draw
    void runFit();
    void saveResults();
    void saveConfig();
    void loadConfigDialog();
    void loadDataDialog();
    void updateEstimateLabels();
    void redraw();
    void setMouseMode(MouseMode m);
    void log(const std::string& line);
    void error(const std::string& message);
    bool dataSettingsChanged(const Config& c) const;

    FitSession session_;
    Config loadedDataSettings_;                  // input / histogram used for the loaded data
    std::string configPath_;
    MouseMode mouseMode_ = kZoom;

    // left
    TRootEmbeddedCanvas* embedded_ = nullptr;
    SpectrumCanvas* canvas_ = nullptr;
    TPad* padLogo_ = nullptr;            // watermark behind the spectrum
    TPad* padSpectrum_ = nullptr;
    TPad* padResiduals_ = nullptr;
    TGRadioButton* modeButtons_[3] = {nullptr, nullptr, nullptr};
    TGCheckButton* logY_ = nullptr;
    TGCheckButton* showComponents_ = nullptr;
    TGLabel* status_ = nullptr;
    std::unique_ptr<TH1D> display_;              // data before a fit
    std::unique_ptr<TH1D> residuals_;

    // right: files and data
    TGLabel* dataLabel_ = nullptr;
    TGLabel* configLabel_ = nullptr;
    TGComboBox* format_ = nullptr;
    TGNumberEntry* column_ = nullptr;
    TGTextEntry* rootObject_ = nullptr;
    TGNumberEntry* binWidth_ = nullptr;
    TGNumberEntry* rangeMin_ = nullptr;
    TGNumberEntry* rangeMax_ = nullptr;

    // right: range and estimates
    TGCheckButton* useFitStart_ = nullptr;
    TGNumberEntry* fitStart_ = nullptr;
    TGNumberEntry* fitStartFraction_ = nullptr;
    TGNumberEntry* fitEnd_ = nullptr;
    TGCheckButton* useBgRange_ = nullptr;
    TGNumberEntry* bgFrom_ = nullptr;
    TGNumberEntry* bgTo_ = nullptr;
    TGLabel* estimates1_ = nullptr;
    TGLabel* estimates2_ = nullptr;

    // right: components
    TGCompositeFrame* resolutionRows_ = nullptr;
    TGCompositeFrame* lifetimeRows_ = nullptr;
    TGLabel* modelLabel_ = nullptr;
    TGNumberEntry* removeIndex_ = nullptr;
    std::vector<ResolutionRow> resRows_;
    std::vector<LifetimeRow> lifeRows_;

    // right: fit settings and results
    TGNumberEntry* iterations_ = nullptr;
    TGNumberEntry* partial_ = nullptr;
    TGTextEntry* firstOptions_ = nullptr;
    TGTextEntry* options_ = nullptr;
    TGCheckButton* profileWidths_ = nullptr;
    TGCheckButton* distribution_ = nullptr;
    TGTextView* report_ = nullptr;
};

} // namespace pals

#endif
