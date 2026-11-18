#include "PalsGui.h"

#include "Spectrum.h"
#include "TextReport.h"

#include <Buttons.h>
#include <TApplication.h>
#include <TBox.h>
#include <TColor.h>
#include <TF1.h>
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGFileDialog.h>
#include <TGLabel.h>
#include <TGMsgBox.h>
#include <TGNumberEntry.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGView.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TRootEmbeddedCanvas.h>
#include <TGCanvas.h>
#include <TGFont.h>
#include <TGIcon.h>
#include <TGPicture.h>
#include <TImage.h>
#include <TGGC.h>
#include <TList.h>
#include <TSystem.h>
#include <TVirtualX.h>
#include <WidgetMessageTypes.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>

namespace pals {

namespace {

enum WidgetId {
    kLoadData = 100, kLoadConfig, kSaveConfig, kReloadData, kEstimate, kAddResolution, kRemoveResolution,
    kAddLifetime, kRemoveLifetime, kFit, kSaveResults, kUnzoom,
    kModeZoom = 200, kModeRange, kModeBackground,
    kLogY = 300, kShowComponents, kUseFitStart, kUseBgRange,
    kFitStart = 400, kFitStartFraction, kFitEnd, kBgFrom, kBgTo
};

// configuration shown before the user loads one
const char* kDefaultConfig = R"({
  "description": "default configuration of the GUI",
  "input": {"format": "times", "column": 0},
  "histogram": {"bin_width": 0.02, "range": [0, 120]},
  "estimation": {"fit_end": 40},
  "resolution": [{"sigma": "auto", "intensity": 1.0}],
  "lifetimes": [
    {"tau": 0.16, "intensity": 0.3, "type": "nf"},
    {"tau": 0.4, "intensity": 0.5, "type": "nf"},
    {"tau": 2.0, "intensity": 0.2, "type": "nf"}
  ]
})";

TGLayoutHints* hints(ULong_t h = kLHintsLeft | kLHintsCenterY, Int_t l = 2, Int_t r = 2, Int_t t = 1, Int_t b = 1)
{
    return new TGLayoutHints(h, l, r, t, b);
}

int modeIndex(ParMode m) { return m == ParMode::Free ? 0 : m == ParMode::Fixed ? 1 : 2; }
ParMode modeFromIndex(int i) { return i == 1 ? ParMode::Fixed : i == 2 ? ParMode::Partial : ParMode::Free; }

TGComboBox* modeCombo(TGCompositeFrame* parent, ParMode m)
{
    auto* cb = new TGComboBox(parent);
    cb->AddEntry("free", 0);
    cb->AddEntry("fixed", 1);
    cb->AddEntry("partial", 2);
    cb->Resize(68, 20);
    cb->Select(modeIndex(m), kFALSE);
    parent->AddFrame(cb, hints());
    return cb;
}

TGNumberEntry* number(TGCompositeFrame* parent, double value, int id = -1, int digits = 6, bool integer = false)
{
    auto* e = new TGNumberEntry(parent, value, digits, id,
                                integer ? TGNumberFormat::kNESInteger : TGNumberFormat::kNESReal,
                                TGNumberFormat::kNEAAnyNumber);
    parent->AddFrame(e, hints());
    return e;
}

TGLabel* label(TGCompositeFrame* parent, const char* text, UInt_t width = 0)
{
    auto* l = new TGLabel(parent, text);
    l->SetTextJustify(kTextLeft);
    if (width) {
        l->SetWidth(width);
        l->ChangeOptions(l->GetOptions() | kFixedWidth);
    }
    parent->AddFrame(l, hints());
    return l;
}

RoundButton* button(TGCompositeFrame* parent, const char* text, int id, const TGWindow* receiver)
{
    auto* b = new RoundButton(parent, text, id);
    b->Associate(receiver);
    parent->AddFrame(b, hints());
    return b;
}

// Path of a file in resources/ (logo images): $PALS_RESOURCES, the source tree (CMake),
// next to the executable, or the working directory. Empty if not found.
std::string resourcePath(const std::string& name)
{
    std::vector<std::string> dirs;
    if (const char* env = gSystem->Getenv("PALS_RESOURCES")) dirs.emplace_back(env);
#ifdef PALS_RESOURCE_DIR
    dirs.emplace_back(PALS_RESOURCE_DIR);
#endif
    if (gApplication && gApplication->Argc() > 0) {
        const std::string exeDir = gSystem->GetDirName(gApplication->Argv(0)).Data();
        dirs.push_back(exeDir + "/resources");
        dirs.push_back(exeDir + "/../resources");
    }
    dirs.emplace_back("resources");
    dirs.emplace_back("../resources");
    for (const auto& d : dirs) {
        const std::string path = d + "/" + name;
        if (!gSystem->AccessPathName(path.c_str())) return path;   // kFALSE = the file exists
    }
    return {};
}

// Picture frame with a logo image, or nullptr if the image is not found.
TGIcon* logoIcon(TGCompositeFrame* parent, const std::string& name)
{
    const std::string path = resourcePath(name);
    if (path.empty()) return nullptr;
    const TGPicture* pic = gClient->GetPicture(path.c_str());
    if (!pic) return nullptr;
    return new TGIcon(parent, pic, pic->GetWidth(), pic->GetHeight());
}

std::string fmt(const char* f, double a, double b = 0.0, double c = 0.0)
{
    char buf[256];
    std::snprintf(buf, sizeof buf, f, a, b, c);
    return buf;
}

} // namespace

// ----------------------------------------------------------------------------- canvas
// ----------------------------------------------------------------------------- look
Pixel_t GuiTheme::panel() { return TColor::RGB2Pixel(0x87, 0xce, 0xeb); }
Pixel_t GuiTheme::text() { return TColor::RGB2Pixel(0x1f, 0x35, 0x52); }
Pixel_t GuiTheme::button() { return TColor::RGB2Pixel(0xf0, 0xf8, 0xfc); }
Pixel_t GuiTheme::buttonHover() { return TColor::RGB2Pixel(0xff, 0xff, 0xff); }
Pixel_t GuiTheme::buttonDown() { return TColor::RGB2Pixel(0xc4, 0xe4, 0xf3); }
Pixel_t GuiTheme::buttonDisabled() { return TColor::RGB2Pixel(0xd3, 0xdc, 0xe1); }
Pixel_t GuiTheme::buttonBorder() { return TColor::RGB2Pixel(0x2e, 0x4e, 0x73); }

FontStruct_t GuiTheme::font(bool bold)
{
    static FontStruct_t cache[2] = {0, 0};
    static bool searched[2] = {false, false};
    const int k = bold ? 1 : 0;
    if (!searched[k]) {
        searched[k] = true;
        const char* weight = bold ? "bold" : "medium";
        for (const char* family : {"helvetica", "nimbus sans", "liberation sans", "dejavu sans", "arial"}) {
            const std::string name = std::string("-*-") + family + "-" + weight + "-r-*-*-12-*-*-*-*-*-*-*";
            if (TGFont* f = gClient->GetFont(name.c_str(), kFALSE)) {
                cache[k] = f->GetFontStruct();
                break;
            }
        }
    }
    return cache[k];
}

RoundButton::RoundButton(const TGWindow* parent, const char* text, Int_t id)
    : TGTextButton(parent, text, id), gc_(new TGGC(TGFrame::GetBckgndGC())), surround_(GetBackground())
{
}

RoundButton::~RoundButton() { delete gc_; }

void RoundButton::ChangeBackground(Pixel_t back)
{
    surround_ = back;
    TGTextButton::ChangeBackground(back);
}

TGDimension RoundButton::GetDefaultSize() const
{
    const TGDimension d = TGTextButton::GetDefaultSize();
    return TGDimension(d.fWidth + GuiTheme::kButtonExtraWidth, d.fHeight + GuiTheme::kButtonExtraHeight);
}

Bool_t RoundButton::HandleCrossing(Event_t* event)
{
    const Bool_t r = TGTextButton::HandleCrossing(event);
    const bool inside = event->fType == kEnterNotify;
    if (inside != hover_) {
        hover_ = inside;
        gClient->NeedRedraw(this);
    }
    return r;
}

void RoundButton::DoRedraw()
{
    const int w = static_cast<int>(fWidth), h = static_cast<int>(fHeight);
    gVirtualX->SetWindowBackground(fId, surround_);
    gVirtualX->ClearArea(fId, 0, 0, fWidth, fHeight);
    if (w < 4 || h < 4) return;

    // rounded rectangle as a polygon (6 segments per corner), inset by d pixels
    auto rounded = [&](int d) {
        const int r = std::max(1, std::min(GuiTheme::kCornerRadius - d, std::min(w, h) / 2 - 1 - d));
        const int x0 = d, y0 = d, x1 = w - 1 - d, y1 = h - 1 - d;
        const int cx[4] = {x1 - r, x0 + r, x0 + r, x1 - r}, cy[4] = {y0 + r, y0 + r, y1 - r, y1 - r};
        std::vector<Point_t> pts;
        for (int c = 0; c < 4; ++c)
            for (int k = 0; k <= 6; ++k) {
                const double a = (c * 90.0 + k * 15.0) * M_PI / 180.0;
                Point_t q;
                q.fX = static_cast<Short_t>(std::lround(cx[c] + r * std::cos(a)));
                q.fY = static_cast<Short_t>(std::lround(cy[c] - r * std::sin(a)));
                pts.push_back(q);
            }
        return pts;
    };

    const bool down = fState == kButtonDown || fState == kButtonEngaged;
    Pixel_t face = GuiTheme::button();
    if (fState == kButtonDisabled) face = GuiTheme::buttonDisabled();
    else if (down) face = GuiTheme::buttonDown();
    else if (hover_) face = GuiTheme::buttonHover();

    const int border = std::min(GuiTheme::kBorderWidth, std::min(w, h) / 4);
    auto outer = rounded(0);
    gc_->SetForeground(GuiTheme::buttonBorder());
    gVirtualX->FillPolygon(fId, gc_->GetGC(), outer.data(), static_cast<Int_t>(outer.size()));
    auto inner = rounded(border);
    gc_->SetForeground(face);
    gVirtualX->FillPolygon(fId, gc_->GetGC(), inner.data(), static_cast<Int_t>(inner.size()));

    // label (position as in TGTextButton)
    int x, y;
    if (fTMode & kTextLeft) x = fMLeft + 4;
    else if (fTMode & kTextRight) x = w - static_cast<int>(fTWidth) - fMRight - 4;
    else x = (w - static_cast<int>(fTWidth) + fMLeft - fMRight) / 2;
    if (fTMode & kTextTop) y = fMTop + 3;
    else if (fTMode & kTextBottom) y = h - static_cast<int>(fTHeight) - fMBottom - 3;
    else y = (h - static_cast<int>(fTHeight) + fMTop - fMBottom) / 2;
    if (down) { ++x; ++y; }
    fTLayout->DrawText(fId, fState == kButtonDisabled ? TGFrame::GetShadowGC()() : fNormGC, x, y, 0, -1);
}

namespace {

// Paints the control area; entry fields, lists, the report and the canvas keep their colours.
void paintFrame(TGFrame* f, Pixel_t color)
{
    if (!f) return;
    if (dynamic_cast<TGTextEntry*>(f) || dynamic_cast<TGNumberEntry*>(f) || dynamic_cast<TGComboBox*>(f) ||
        dynamic_cast<TGCanvas*>(f) || dynamic_cast<TGView*>(f))
        return;
    f->ChangeBackground(color);
    if (auto* c = dynamic_cast<TGCompositeFrame*>(f)) {
        if (!c->GetList()) return;
        TIter next(c->GetList());
        while (auto* el = static_cast<TGFrameElement*>(next())) paintFrame(el->fFrame, color);
    }
}

// Navy-grey sans-serif text on labels, buttons, group titles and entry fields
// (the report keeps its fixed-width font, which the tables need).
void styleText(TGFrame* f)
{
    if (!f || dynamic_cast<TGComboBox*>(f) || dynamic_cast<TGCanvas*>(f) || dynamic_cast<TGView*>(f)) return;
    const Pixel_t color = GuiTheme::text();
    const FontStruct_t normal = GuiTheme::font(false), bold = GuiTheme::font(true);
    if (auto* l = dynamic_cast<TGLabel*>(f)) {
        if (normal) l->SetTextFont(normal);
        l->SetTextColor(color);
    } else if (auto* b = dynamic_cast<TGTextButton*>(f)) {
        // push buttons bold, check / radio buttons normal
        const FontStruct_t font = dynamic_cast<RoundButton*>(f) && bold ? bold : normal;
        if (font) b->SetFont(font);
        b->SetTextColor(color);
    } else if (auto* e = dynamic_cast<TGTextEntry*>(f)) {
        if (normal) e->SetFont(normal);
        e->SetTextColor(color);
    } else if (auto* g = dynamic_cast<TGGroupFrame*>(f)) {
        if (bold) g->SetTextFont(bold);
        else if (normal) g->SetTextFont(normal);
        g->SetTextColor(color);
    }
    if (auto* c = dynamic_cast<TGCompositeFrame*>(f)) {
        if (!c->GetList()) return;
        TIter next(c->GetList());
        while (auto* el = static_cast<TGFrameElement*>(next())) styleText(el->fFrame);
    }
}

void applyThemeTo(TGFrame* f)
{
    paintFrame(f, GuiTheme::panel());
    styleText(f);
}

} // namespace

void PalsGui::applyTheme()
{
    applyThemeTo(this);
    gClient->NeedRedraw(this, kTRUE);
}

SpectrumCanvas::SpectrumCanvas(const char* name, Int_t ww, Int_t wh, Int_t winid, PalsGui* gui)
    : TCanvas(name, ww, wh, winid), gui_(gui)
{
}

void SpectrumCanvas::HandleInput(EEventType event, Int_t px, Int_t py)
{
    TPad* pad = gui_ ? gui_->spectrumPad() : nullptr;
    if (!pad || gui_->mouseMode() == PalsGui::kZoom) {
        TCanvas::HandleInput(event, px, py);
        if (pad && event == kMouseMotion) gui_->cursorMoved(pad->AbsPixeltoX(px));
        return;
    }
    const double x = pad->AbsPixeltoX(px);
    switch (event) {
        case kButton1Down:
            selecting_ = true;
            x0_ = x;
            break;
        case kButton1Motion:
            if (selecting_) gui_->selectionMoved(x0_, x);
            break;
        case kButton1Up:
            if (selecting_) {
                selecting_ = false;
                gui_->rangeSelected(x0_, x);
            }
            break;
        case kMouseMotion:
            gui_->cursorMoved(x);
            break;
        default:
            TCanvas::HandleInput(event, px, py);
    }
}

// ----------------------------------------------------------------------------- construction
PalsGui::PalsGui(const TGWindow* parent, UInt_t w, UInt_t h) : TGMainFrame(parent, w, h, kHorizontalFrame)
{
    SetCleanup(kDeepCleanup);
    auto* left = new TGVerticalFrame(this, w > 620 ? w - 620 : 400, h);
    AddFrame(left, hints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
    auto* right = new TGVerticalFrame(this, 610, h, kFixedWidth);
    AddFrame(right, hints(kLHintsRight | kLHintsExpandY, 2, 2, 2, 2));
    buildLeft(left);
    buildRight(right);

    SetWindowName("PALS_avalanche");
    {
        const std::string icon = resourcePath("pav_logo.png");
        if (!icon.empty()) SetIconPixmap(icon.c_str());
    }
    applyTheme();
    MapSubwindows();
    Resize(w, h);
    MapWindow();

    const Config cfg = Config::parse(kDefaultConfig);
    session_.setConfig(cfg);
    loadedDataSettings_ = cfg;
    populate(cfg);
    redraw();
}

PalsGui::~PalsGui() = default;

void PalsGui::CloseWindow()
{
    gApplication->Terminate(0);
}

TGGroupFrame* PalsGui::group(TGCompositeFrame* parent, const char* title)
{
    auto* g = new TGGroupFrame(parent, title);
    parent->AddFrame(g, hints(kLHintsExpandX | kLHintsTop, 2, 2, 2, 2));
    return g;
}

void PalsGui::buildLeft(TGCompositeFrame* parent)
{
    auto* bar = new TGHorizontalFrame(parent);
    parent->AddFrame(bar, hints(kLHintsExpandX | kLHintsTop));
    label(bar, "Mouse:");
    const char* names[3] = {"zoom", "fit range", "background window"};
    for (int i = 0; i < 3; ++i) {
        modeButtons_[i] = new TGRadioButton(bar, names[i], kModeZoom + i);
        modeButtons_[i]->Associate(this);
        bar->AddFrame(modeButtons_[i], hints());
    }
    modeButtons_[0]->SetState(kButtonDown);
    logY_ = new TGCheckButton(bar, "log y", kLogY);
    logY_->Associate(this);
    logY_->SetOn(kTRUE);
    bar->AddFrame(logY_, hints(kLHintsLeft | kLHintsCenterY, 16, 2, 1, 1));
    showComponents_ = new TGCheckButton(bar, "components", kShowComponents);
    showComponents_->Associate(this);
    showComponents_->SetOn(kTRUE);
    bar->AddFrame(showComponents_, hints());
    button(bar, "Unzoom", kUnzoom, this);

    embedded_ = new TRootEmbeddedCanvas(nullptr, parent, 900, 800);
    parent->AddFrame(embedded_, hints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
    canvas_ = new SpectrumCanvas("pals_canvas", 10, 10, embedded_->GetCanvasWindowId(), this);
    embedded_->AdoptCanvas(canvas_);
    canvas_->cd();
    if (GuiTheme::kLogoWatermark) {
        const std::string path = resourcePath("pav_logo_watermark.png");
        if (!path.empty()) {
            if (TImage* img = TImage::Open(path.c_str())) {
                // centred in the plot area of the spectrum pad, drawn below it
                padLogo_ = new TPad("logo", "logo", 0.30, 0.42, 0.78, 0.90);
                padLogo_->SetBorderMode(0);
                padLogo_->SetMargin(0.0, 0.0, 0.0, 0.0);
                padLogo_->Draw();
                padLogo_->cd();
                img->SetConstRatio(kTRUE);
                img->Draw();
                padLogo_->SetEditable(kFALSE);
                canvas_->cd();
            }
        }
    }
    padSpectrum_ = new TPad("spectrum", "spectrum", 0.0, 0.28, 1.0, 1.0);
    if (padLogo_) {                       // transparent, so that the watermark shows through
        padSpectrum_->SetFillStyle(4000);
        padSpectrum_->SetFrameFillStyle(4000);
    }
    padSpectrum_->SetBottomMargin(0.05);
    padSpectrum_->SetLeftMargin(0.09);
    padSpectrum_->SetRightMargin(0.02);
    padResiduals_ = new TPad("residuals", "residuals", 0.0, 0.0, 1.0, 0.28);
    padResiduals_->SetTopMargin(0.02);
    padResiduals_->SetBottomMargin(0.3);
    padResiduals_->SetLeftMargin(0.09);
    padResiduals_->SetRightMargin(0.02);
    padSpectrum_->Draw();
    padResiduals_->Draw();

    auto* bottom = new TGHorizontalFrame(parent);
    parent->AddFrame(bottom, hints(kLHintsExpandX | kLHintsBottom, 0, 0, 0, 0));
    if (GuiTheme::kLogoBottomLeft)
        if (TGIcon* icon = logoIcon(bottom, "pav_logo_mark.png"))
            bottom->AddFrame(icon, hints(kLHintsLeft | kLHintsCenterY, 4, 6, 2, 2));
    status_ = new TGLabel(bottom, "Load data and a configuration (or edit the settings on the right).");
    status_->SetTextJustify(kTextLeft);
    bottom->AddFrame(status_, hints(kLHintsExpandX | kLHintsCenterY, 4, 4, 2, 2));
}

void PalsGui::buildRight(TGCompositeFrame* parent)
{
    // --- files
    auto* files = group(parent, "Files");
    auto* row = new TGHorizontalFrame(files);
    files->AddFrame(row, hints(kLHintsExpandX));
    button(row, "Load data...", kLoadData, this);
    button(row, "Load config...", kLoadConfig, this);
    button(row, "Save config...", kSaveConfig, this);
    dataLabel_ = label(files, "data: -", 580);
    configLabel_ = label(files, "config: (default)", 580);

    // --- data
    auto* data = group(parent, "Data");
    row = new TGHorizontalFrame(data);
    data->AddFrame(row, hints(kLHintsExpandX));
    label(row, "format");
    format_ = new TGComboBox(row);
    format_->AddEntry("times", 0);
    format_->AddEntry("histogram", 1);
    format_->AddEntry("ROOT", 2);
    format_->Resize(90, 20);
    row->AddFrame(format_, hints());
    label(row, "column");
    column_ = number(row, 0, -1, 3, true);
    label(row, "ROOT object");
    rootObject_ = new TGTextEntry(row, "");
    rootObject_->Resize(170, 20);
    row->AddFrame(rootObject_, hints());
    row = new TGHorizontalFrame(data);
    data->AddFrame(row, hints(kLHintsExpandX));
    label(row, "bin width [ns]");
    binWidth_ = number(row, 0.02);
    label(row, "range [ns]");
    rangeMin_ = number(row, 0);
    rangeMax_ = number(row, 120);
    button(row, "Reload", kReloadData, this);

    // --- fit range and estimates
    auto* range = group(parent, "Fit range and estimates (Enter or mouse selection updates)");
    row = new TGHorizontalFrame(range);
    range->AddFrame(row, hints(kLHintsExpandX));
    useFitStart_ = new TGCheckButton(row, "start [ns]", kUseFitStart);
    useFitStart_->Associate(this);
    row->AddFrame(useFitStart_, hints());
    fitStart_ = number(row, 0, kFitStart);
    fitStart_->Associate(this);
    label(row, "or fraction of peak");
    fitStartFraction_ = number(row, 0.3, kFitStartFraction, 4);
    fitStartFraction_->Associate(this);
    label(row, "end [ns]");
    fitEnd_ = number(row, 40, kFitEnd);
    fitEnd_->Associate(this);
    row = new TGHorizontalFrame(range);
    range->AddFrame(row, hints(kLHintsExpandX));
    useBgRange_ = new TGCheckButton(row, "background window [ns]", kUseBgRange);
    useBgRange_->Associate(this);
    row->AddFrame(useBgRange_, hints());
    bgFrom_ = number(row, 0, kBgFrom);
    bgFrom_->Associate(this);
    bgTo_ = number(row, 0, kBgTo);
    bgTo_->Associate(this);
    button(row, "Estimate", kEstimate, this);
    estimates1_ = label(range, "-", 580);
    estimates2_ = label(range, "-", 580);

    // --- resolution
    auto* res = group(parent, "Resolution (Gaussians)");
    row = new TGHorizontalFrame(res);
    res->AddFrame(row, hints(kLHintsExpandX));
    label(row, "sigma [ns]", 80);
    label(row, "auto", 44);
    label(row, "intensity", 80);
    label(row, "sigma", 72);
    label(row, "intensity", 72);
    resolutionRows_ = new TGVerticalFrame(res);
    resolutionRows_->SetCleanup(kDeepCleanup);
    res->AddFrame(resolutionRows_, hints(kLHintsExpandX));
    row = new TGHorizontalFrame(res);
    res->AddFrame(row, hints(kLHintsExpandX));
    button(row, "+ add", kAddResolution, this);
    button(row, "- remove last", kRemoveResolution, this);

    // --- lifetimes
    auto* life = group(parent, "Lifetime components");
    row = new TGHorizontalFrame(life);
    life->AddFrame(row, hints(kLHintsExpandX));
    label(row, "#", 18);
    label(row, "tau [ns]", 80);
    label(row, "intensity", 80);
    label(row, "tau", 72);
    label(row, "intensity", 72);
    label(row, "p-Ps", 44);
    label(row, "width", 132);
    label(row, "width", 60);
    lifetimeRows_ = new TGVerticalFrame(life);
    lifetimeRows_->SetCleanup(kDeepCleanup);
    life->AddFrame(lifetimeRows_, hints(kLHintsExpandX));
    row = new TGHorizontalFrame(life);
    life->AddFrame(row, hints(kLHintsExpandX));
    button(row, "+ add", kAddLifetime, this);
    button(row, "- remove #", kRemoveLifetime, this);
    removeIndex_ = number(row, 1, -1, 3, true);
    modelLabel_ = label(life, "p-Ps: intensity = ratio to o-Ps; width = sd of ln(tau)", 580);

    // --- fit settings
    auto* fitGroup = group(parent, "Fit");
    auto* fitBox = new TGHorizontalFrame(fitGroup);
    fitGroup->AddFrame(fitBox, hints(kLHintsExpandX, 0, 0, 0, 0));
    auto* fit = new TGVerticalFrame(fitBox);
    fitBox->AddFrame(fit, hints(kLHintsLeft | kLHintsExpandX | kLHintsCenterY, 0, 0, 0, 0));
    if (GuiTheme::kLogoInFitPanel)
        if (TGIcon* icon = logoIcon(fitBox, "pav_logo_card.png"))
            fitBox->AddFrame(icon, hints(kLHintsRight | kLHintsCenterY, 6, 2, 0, 0));
    row = new TGHorizontalFrame(fit);
    fit->AddFrame(row, hints(kLHintsExpandX));
    label(row, "iterations");
    iterations_ = number(row, 3, -1, 3, true);
    label(row, "partial %");
    partial_ = number(row, 5, -1, 4);
    label(row, "options: first");
    firstOptions_ = new TGTextEntry(row, "RMW");
    firstOptions_->Resize(50, 20);
    row->AddFrame(firstOptions_, hints());
    label(row, "then");
    options_ = new TGTextEntry(row, "RM");
    options_->Resize(50, 20);
    row->AddFrame(options_, hints());
    row = new TGHorizontalFrame(fit);
    fit->AddFrame(row, hints(kLHintsExpandX));
    profileWidths_ = new TGCheckButton(row, "profile intervals of widths");
    row->AddFrame(profileWidths_, hints());
    distribution_ = new TGCheckButton(row, "distribution analysis (NNLS bagged + MaxEnt)");
    row->AddFrame(distribution_, hints());

    // --- actions and report
    row = new TGHorizontalFrame(parent);
    parent->AddFrame(row, hints(kLHintsExpandX, 2, 2, 4, 4));
    auto* fitButton = new RoundButton(row, "  Fit  ", kFit);
    fitButton->Associate(this);
    row->AddFrame(fitButton, hints(kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2));
    auto* saveButton = new RoundButton(row, "Save results...", kSaveResults);
    saveButton->Associate(this);
    row->AddFrame(saveButton, hints(kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2));

    report_ = new TGTextView(parent, 600, 200);
    parent->AddFrame(report_, hints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
}

// ----------------------------------------------------------------------------- rows
void PalsGui::rebuildResolutionRows(const std::vector<ResolutionComponent>& rows)
{
    resolutionRows_->Cleanup();
    resRows_.clear();
    for (const auto& r : rows) {
        ResolutionRow w{};
        w.frame = new TGHorizontalFrame(resolutionRows_);
        resolutionRows_->AddFrame(w.frame, hints(kLHintsExpandX, 0, 0, 0, 0));
        w.sigma = number(w.frame, r.sigmaAuto ? 0.1 : r.sigma);
        w.sigmaAuto = new TGCheckButton(w.frame, "");
        w.sigmaAuto->SetOn(r.sigmaAuto);
        w.frame->AddFrame(w.sigmaAuto, hints(kLHintsLeft | kLHintsCenterY, 14, 14, 1, 1));
        w.intensity = number(w.frame, r.intensity);
        w.sigmaMode = modeCombo(w.frame, r.sigmaMode);
        w.intensityMode = modeCombo(w.frame, r.intensityMode);
        resRows_.push_back(w);
    }
    applyThemeTo(resolutionRows_);
    resolutionRows_->MapSubwindows();
    resolutionRows_->Layout();
    MapSubwindows();
    Layout();
}

void PalsGui::rebuildLifetimeRows(const std::vector<LifetimeComponent>& rows)
{
    lifetimeRows_->Cleanup();
    lifeRows_.clear();
    int k = 0;
    for (const auto& l : rows) {
        LifetimeRow w{};
        w.frame = new TGHorizontalFrame(lifetimeRows_);
        lifetimeRows_->AddFrame(w.frame, hints(kLHintsExpandX, 0, 0, 0, 0));
        label(w.frame, std::to_string(++k).c_str(), 14);
        w.tau = number(w.frame, l.tau);
        w.intensity = number(w.frame, l.intensity);
        w.tauMode = modeCombo(w.frame, l.tauMode);
        w.intensityMode = modeCombo(w.frame, l.pPs ? ParMode::Fixed : l.intensityMode);
        w.pPs = new TGCheckButton(w.frame, "");
        w.pPs->SetOn(l.pPs);
        w.frame->AddFrame(w.pPs, hints(kLHintsLeft | kLHintsCenterY, 14, 14, 1, 1));
        w.hasWidth = new TGCheckButton(w.frame, "");
        w.hasWidth->SetOn(l.hasWidth);
        w.frame->AddFrame(w.hasWidth, hints());
        w.width = number(w.frame, l.hasWidth ? l.width : 0.05, -1, 5);
        w.widthMode = modeCombo(w.frame, l.widthMode);
        lifeRows_.push_back(w);
    }
    applyThemeTo(lifetimeRows_);
    lifetimeRows_->MapSubwindows();
    lifetimeRows_->Layout();
    MapSubwindows();
    Layout();
}

std::vector<ResolutionComponent> PalsGui::currentResolution() const
{
    std::vector<ResolutionComponent> out;
    for (const auto& w : resRows_) {
        ResolutionComponent r;
        r.sigmaAuto = w.sigmaAuto->IsOn();
        r.sigma = w.sigma->GetNumber();
        r.intensity = w.intensity->GetNumber();
        r.sigmaMode = modeFromIndex(w.sigmaMode->GetSelected());
        r.intensityMode = modeFromIndex(w.intensityMode->GetSelected());
        out.push_back(r);
    }
    return out;
}

std::vector<LifetimeComponent> PalsGui::currentLifetimes() const
{
    std::vector<LifetimeComponent> out;
    for (const auto& w : lifeRows_) {
        LifetimeComponent l;
        l.tau = w.tau->GetNumber();
        l.intensity = w.intensity->GetNumber();
        l.tauMode = modeFromIndex(w.tauMode->GetSelected());
        l.pPs = w.pPs->IsOn();
        l.intensityMode = l.pPs ? ParMode::Fixed : modeFromIndex(w.intensityMode->GetSelected());
        l.hasWidth = !l.pPs && w.hasWidth->IsOn();
        l.width = w.width->GetNumber();
        l.widthMode = modeFromIndex(w.widthMode->GetSelected());
        out.push_back(l);
    }
    return out;
}

// ----------------------------------------------------------------------------- config <-> widgets
void PalsGui::populate(const Config& c)
{
    format_->Select(c.input.format == InputFormat::Times ? 0 : c.input.format == InputFormat::Histogram ? 1 : 2, kFALSE);
    column_->SetNumber(c.input.column);
    rootObject_->SetText(c.input.rootObject.c_str(), kFALSE);
    binWidth_->SetNumber(c.histogram.binWidth);
    rangeMin_->SetNumber(c.histogram.min);
    rangeMax_->SetNumber(c.histogram.max);

    const auto& e = c.estimation;
    useFitStart_->SetOn(e.useFitStart);
    fitStart_->SetNumber(e.useFitStart ? e.fitStart : 0.0);
    fitStartFraction_->SetNumber(e.fitStartFraction);
    fitEnd_->SetNumber(e.fitEnd);
    useBgRange_->SetOn(e.background.useRange);
    bgFrom_->SetNumber(e.background.from);
    bgTo_->SetNumber(e.background.to);

    rebuildResolutionRows(c.resolution);
    if (c.model.type == ModelType::Formula) {
        rebuildLifetimeRows({});
        modelLabel_->SetText(("Formula model (edit it in the config file): " + c.model.formula).c_str());
    } else {
        rebuildLifetimeRows(c.lifetimes);
        modelLabel_->SetText("p-Ps: intensity = ratio to o-Ps; width = sd of ln(tau)");
    }

    iterations_->SetNumber(c.fit.iterations);
    partial_->SetNumber(100.0 * c.fit.partialVariation);
    firstOptions_->SetText(c.fit.firstPassOptions.c_str(), kFALSE);
    options_->SetText(c.fit.options.c_str(), kFALSE);
    profileWidths_->SetOn(c.fit.profileWidths);
    distribution_->SetOn(c.distribution.enabled);
    Layout();
}

Config PalsGui::collect() const
{
    Config c = session_.config();
    const int f = format_->GetSelected();
    c.input.format = f == 1 ? InputFormat::Histogram : f == 2 ? InputFormat::Root : InputFormat::Times;
    c.input.column = static_cast<int>(std::lround(column_->GetNumber()));
    c.input.rootObject = rootObject_->GetText();
    c.histogram.binWidth = binWidth_->GetNumber();
    c.histogram.min = rangeMin_->GetNumber();
    c.histogram.max = rangeMax_->GetNumber();

    auto& e = c.estimation;
    e.useFitStart = useFitStart_->IsOn();
    e.fitStart = fitStart_->GetNumber();
    e.fitStartFraction = fitStartFraction_->GetNumber();
    e.fitEnd = fitEnd_->GetNumber();
    e.background.useRange = useBgRange_->IsOn();
    e.background.from = bgFrom_->GetNumber();
    e.background.to = bgTo_->GetNumber();

    if (c.model.type != ModelType::Formula || c.model.convolve) c.resolution = currentResolution();
    if (c.model.type == ModelType::Lifetimes) c.lifetimes = currentLifetimes();

    c.fit.iterations = static_cast<int>(std::lround(iterations_->GetNumber()));
    c.fit.partialVariation = partial_->GetNumber() / 100.0;
    c.fit.firstPassOptions = firstOptions_->GetText();
    c.fit.options = options_->GetText();
    c.fit.profileWidths = profileWidths_->IsOn();
    c.distribution.enabled = distribution_->IsOn();

    // validation: the same checks as for a config file
    Config v = Config::parse(c.toJson(), configPath_);
    v.sourceText = v.toJson();
    return v;
}

bool PalsGui::dataSettingsChanged(const Config& c) const
{
    const auto& a = c.input;
    const auto& b = loadedDataSettings_.input;
    const auto& ha = c.histogram;
    const auto& hb = loadedDataSettings_.histogram;
    return a.format != b.format || a.column != b.column || a.rootObject != b.rootObject ||
           ha.binWidth != hb.binWidth || ha.min != hb.min || ha.max != hb.max;
}

// ----------------------------------------------------------------------------- actions
void PalsGui::openConfig(const std::string& path)
{
    const Config c = Config::load(path);
    configPath_ = path;
    session_.setConfig(c, path);
    populate(c);
    configLabel_->SetText(("config: " + path).c_str());
    if (session_.hasData()) {
        if (dataSettingsChanged(c)) {
            session_.loadData(session_.dataPath());
            loadedDataSettings_ = c;
            display_.reset();
        }
        applyAndEstimate();
    } else {
        loadedDataSettings_ = c;
        redraw();
    }
    log("Config loaded: " + path);
}

void PalsGui::openData(const std::string& path)
{
    const Config c = collect();
    session_.setConfig(c, configPath_);
    session_.loadData(path);
    loadedDataSettings_ = c;
    display_.reset();
    residuals_.reset();
    dataLabel_->SetText(("data: " + path).c_str());
    log("Data loaded: " + path);
    applyAndEstimate();
}

void PalsGui::applyAndEstimate()
{
    const Config c = collect();
    const bool reload = session_.hasData() && dataSettingsChanged(c);
    session_.setConfig(c, configPath_);
    residuals_.reset();
    if (!session_.hasData()) {
        redraw();
        return;
    }
    if (reload) {
        session_.loadData(session_.dataPath());
        loadedDataSettings_ = c;
        display_.reset();
        log("Data re-read with the new data settings");
    }
    session_.estimate();
    updateEstimateLabels();
    redraw();
}

void PalsGui::runFit()
{
    if (!session_.hasData()) {
        error("Load data first.");
        return;
    }
    applyAndEstimate();
    report_->Clear();
    log("Fitting... (ROOT fit output in the terminal)");
    session_.fit([this](const StageSummary& s) {
        log(fmt("stage %.0f: status %.0f, chi2/ndf = %.4f", s.stage, s.status, s.ndf > 0 ? s.chi2 / s.ndf : 0.0));
        gSystem->ProcessEvents();
    });

    // residuals (data - fit)/sqrt(data) in the fit range
    const Spectrum& sp = session_.spectrum();
    const Estimates& est = session_.estimates();
    residuals_ = std::make_unique<TH1D>("residuals", ";time [ns];(data-fit)/#sqrt{data}", static_cast<int>(sp.size()),
                                        sp.min, sp.max());
    residuals_->SetDirectory(nullptr);
    residuals_->SetStats(kFALSE);
    TF1* f = session_.function();
    for (std::size_t i = est.fitFrom; i <= est.fitTo; ++i)
        if (sp.counts[i] > 0)
            residuals_->SetBinContent(static_cast<int>(i) + 1, (sp.counts[i] - f->Eval(sp.center(i))) / std::sqrt(sp.counts[i]));

    std::ostringstream os;
    session_.writeReport(os);
    report_->LoadBuffer(os.str().c_str());
    updateEstimateLabels();
    redraw();
    const auto& o = session_.outcome();
    status_->SetText(fmt("Fit done: chi2/ndf = %.4f, status %.0f", o.ndf > 0 ? o.chi2 / o.ndf : 0.0, o.status).c_str());
}

void PalsGui::saveResults()
{
    if (!session_.hasFit()) {
        error("Run a fit first.");
        return;
    }
    static const char* types[] = {"Results (txt + root)", "*.txt", "All files", "*", nullptr, nullptr};
    TGFileInfo fi;
    fi.fFileTypes = types;
    fi.fFilename = StrDup(session_.defaultOutputBase().c_str());
    new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fi);
    if (!fi.fFilename) return;
    std::string base = fi.fFilename;
    for (const char* ext : {".txt", ".root"})
        if (base.size() > std::string(ext).size() && base.compare(base.size() - std::string(ext).size(), std::string::npos, ext) == 0)
            base.erase(base.size() - std::string(ext).size());
    session_.save(base);
    log("Results saved: " + base + ".txt / .root");
}

void PalsGui::saveConfig()
{
    const Config c = collect();
    static const char* types[] = {"JSON config", "*.json", "All files", "*", nullptr, nullptr};
    TGFileInfo fi;
    fi.fFileTypes = types;
    fi.fFilename = StrDup(configPath_.empty() ? "config.json" : configPath_.c_str());
    new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fi);
    if (!fi.fFilename) return;
    std::string path = fi.fFilename;
    if (path.size() < 5 || path.compare(path.size() - 5, 5, ".json") != 0) path += ".json";
    c.save(path);
    configPath_ = path;          // the session keeps its estimates and fit
    configLabel_->SetText(("config: " + path).c_str());
    log("Config saved: " + path);
}

void PalsGui::loadConfigDialog()
{
    static const char* types[] = {"JSON config", "*.json", "All files", "*", nullptr, nullptr};
    TGFileInfo fi;
    fi.fFileTypes = types;
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (fi.fFilename) openConfig(fi.fFilename);
}

void PalsGui::loadDataDialog()
{
    static const char* types[] = {"All files", "*", "Text", "*.dat", "ROOT", "*.root", nullptr, nullptr};
    TGFileInfo fi;
    fi.fFileTypes = types;
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (fi.fFilename) openData(fi.fFilename);
}

void PalsGui::updateEstimateLabels()
{
    if (!session_.hasEstimates()) {
        estimates1_->SetText("-");
        estimates2_->SetText("-");
        return;
    }
    const Spectrum& s = session_.spectrum();
    const Estimates& e = session_.estimates();
    char buf[256];
    std::snprintf(buf, sizeof buf, "fit %.3f - %.3f ns,  peak %.3f ns,  sigma (left HWHM) %.3f ns", s.center(e.fitFrom),
                  s.center(e.fitTo), e.peakPosition, e.sigma);
    estimates1_->SetText(buf);
    std::snprintf(buf, sizeof buf, "background %.3f +- %.3f counts/bin (%.2f - %.2f ns),  area %.4g counts",
                  e.background, e.backgroundError, s.center(e.backgroundFrom), s.center(e.backgroundTo), e.area);
    estimates2_->SetText(buf);
}

// ----------------------------------------------------------------------------- drawing
void PalsGui::redraw()
{
    padSpectrum_->cd();
    padSpectrum_->Clear();
    padResiduals_->cd();
    padResiduals_->Clear();
    if (!session_.hasData()) {
        canvas_->cd();
        canvas_->Modified();
        canvas_->Update();
        return;
    }
    const Spectrum& s = session_.spectrum();
    if (!display_) {
        display_ = std::make_unique<TH1D>("pals_data", ";time [ns];counts", static_cast<int>(s.size()), s.min, s.max());
        display_->SetDirectory(nullptr);
        display_->SetStats(kFALSE);
        for (std::size_t i = 0; i < s.size(); ++i) display_->SetBinContent(static_cast<int>(i) + 1, s.counts[i]);
        if (session_.hasEstimates()) {
            const Estimates& e = session_.estimates();
            display_->GetXaxis()->SetRangeUser(s.center(e.firstBin) - 1.0, s.center(e.lastBin));
        }
        display_->SetLineColor(kBlack);
    }

    padSpectrum_->cd();
    padSpectrum_->SetLogy(logY_->IsOn() ? 1 : 0);
    const bool logy = logY_->IsOn();
    const double ymax = display_->GetMaximum() * (logy ? 3.0 : 1.1);
    const double ymin = logy ? 0.5 : 0.0;
    display_->SetMinimum(ymin);
    display_->SetMaximum(ymax);
    display_->GetXaxis()->SetLabelSize(0.035);
    display_->Draw("HIST");

    auto* legend = new TLegend(0.72, 0.72, 0.97, 0.95);
    legend->SetBit(kCanDelete);
    legend->AddEntry(display_.get(), "data", "l");

    if (session_.hasEstimates()) {
        const Estimates& e = session_.estimates();
        auto* fitBox = new TBox(s.center(e.fitFrom), ymin, s.center(e.fitTo), ymax);
        fitBox->SetFillStyle(3004);
        fitBox->SetFillColor(kAzure - 9);
        fitBox->SetLineColor(kAzure - 9);
        fitBox->SetBit(kCanDelete);
        fitBox->Draw();
        for (double x : {s.center(e.fitFrom), s.center(e.fitTo)}) {
            auto* l = new TLine(x, ymin, x, ymax);
            l->SetLineColor(kAzure + 2);
            l->SetLineWidth(2);
            l->SetBit(kCanDelete);
            l->Draw();
        }
        auto* bgBox = new TBox(s.center(e.backgroundFrom), ymin, s.center(e.backgroundTo), ymax);
        bgBox->SetFillStyle(3005);
        bgBox->SetFillColor(kOrange + 1);
        bgBox->SetLineColor(kOrange + 1);
        bgBox->SetBit(kCanDelete);
        bgBox->Draw();
        const double B = session_.hasFit() ? session_.outcome().params[0] : e.background;
        auto* bl = new TLine(display_->GetXaxis()->GetBinLowEdge(display_->GetXaxis()->GetFirst()), B,
                             display_->GetXaxis()->GetBinUpEdge(display_->GetXaxis()->GetLast()), B);
        bl->SetLineColor(kGray + 2);
        bl->SetLineStyle(3);
        bl->SetBit(kCanDelete);
        bl->Draw();
        legend->AddEntry(fitBox, "fit range", "f");
        legend->AddEntry(bgBox, "background window", "f");
    }
    display_->Draw("HIST SAME AXIS");
    display_->Draw("HIST SAME");

    if (session_.hasFit()) {
        TF1* f = session_.function();
        f->SetLineColor(kRed);
        f->SetLineWidth(2);
        f->Draw("SAME");
        legend->AddEntry(f, "fit", "l");
        if (showComponents_->IsOn()) {
            for (TF1* c : session_.components()) c->Draw("SAME");
            if (!session_.components().empty()) legend->AddEntry(session_.components().front(), "components", "l");
        }
    }
    legend->Draw();

    padResiduals_->cd();
    if (residuals_) {
        residuals_->GetXaxis()->SetRange(display_->GetXaxis()->GetFirst(), display_->GetXaxis()->GetLast());
        residuals_->GetXaxis()->SetLabelSize(0.10);
        residuals_->GetXaxis()->SetTitleSize(0.11);
        residuals_->GetYaxis()->SetLabelSize(0.09);
        residuals_->GetYaxis()->SetTitleSize(0.09);
        residuals_->GetYaxis()->SetTitleOffset(0.4);
        residuals_->SetLineColor(kGray + 2);
        residuals_->Draw("HIST");
        auto* zero = new TLine(display_->GetXaxis()->GetBinLowEdge(display_->GetXaxis()->GetFirst()), 0.0,
                               display_->GetXaxis()->GetBinUpEdge(display_->GetXaxis()->GetLast()), 0.0);
        zero->SetLineColor(kRed);
        zero->SetBit(kCanDelete);
        zero->Draw();
    }
    canvas_->cd();
    canvas_->Modified();
    canvas_->Update();
}

void PalsGui::setMouseMode(MouseMode m)
{
    mouseMode_ = m;
    for (int i = 0; i < 3; ++i) modeButtons_[i]->SetState(i == m ? kButtonDown : kButtonUp);
    const char* help[3] = {"Mouse: ROOT zoom (drag on an axis)", "Mouse: drag over the spectrum to set the fit range",
                           "Mouse: drag over the spectrum to set the background window"};
    status_->SetText(help[m]);
}

void PalsGui::rangeSelected(double x1, double x2)
{
    if (x2 < x1) std::swap(x1, x2);
    if (!session_.hasData() || x2 - x1 < 2.0 * session_.spectrum().binWidth) return;
    try {
        if (mouseMode_ == kFitRange) {
            useFitStart_->SetOn(kTRUE);
            fitStart_->SetNumber(x1);
            fitEnd_->SetNumber(x2);
        } else if (mouseMode_ == kBackground) {
            useBgRange_->SetOn(kTRUE);
            bgFrom_->SetNumber(x1);
            bgTo_->SetNumber(x2);
        }
        applyAndEstimate();
    } catch (const std::exception& e) {
        error(e.what());
    }
}

void PalsGui::selectionMoved(double x1, double x2)
{
    status_->SetText(fmt("selection %.3f - %.3f ns", std::min(x1, x2), std::max(x1, x2)).c_str());
}

void PalsGui::cursorMoved(double x)
{
    if (!display_) return;
    const int bin = display_->FindBin(x);
    status_->SetText(fmt("t = %.3f ns,  counts = %.0f", x, display_->GetBinContent(bin)).c_str());
}

void PalsGui::log(const std::string& line)
{
    report_->AddLine(line.c_str());
    report_->ShowBottom();
}

void PalsGui::error(const std::string& message)
{
    std::cerr << "ERROR: " << message << "\n";
    new TGMsgBox(gClient->GetRoot(), this, "PALS_avalanche", message.c_str(), kMBIconStop, kMBOk);
}

// ----------------------------------------------------------------------------- events
Bool_t PalsGui::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
    try {
        switch (GET_MSG(msg)) {
            case kC_COMMAND:
                switch (GET_SUBMSG(msg)) {
                    case kCM_BUTTON:
                        switch (parm1) {
                            case kLoadData: loadDataDialog(); break;
                            case kLoadConfig: loadConfigDialog(); break;
                            case kSaveConfig: saveConfig(); break;
                            case kReloadData:
                            case kEstimate: applyAndEstimate(); break;
                            case kAddResolution: {
                                auto rows = currentResolution();
                                ResolutionComponent r;
                                r.sigma = rows.empty() ? 0.1 : 2.0 * rows.back().sigma;
                                r.intensity = rows.empty() ? 1.0 : 0.1;
                                rows.push_back(r);
                                rebuildResolutionRows(rows);
                                break;
                            }
                            case kRemoveResolution: {
                                auto rows = currentResolution();
                                if (rows.size() > 1) rows.pop_back();
                                rebuildResolutionRows(rows);
                                break;
                            }
                            case kAddLifetime: {
                                if (session_.config().model.type == ModelType::Formula) break;
                                auto rows = currentLifetimes();
                                LifetimeComponent l;
                                l.tau = rows.empty() ? 0.2 : 2.0 * rows.back().tau;
                                l.intensity = 0.1;
                                rows.push_back(l);
                                rebuildLifetimeRows(rows);
                                break;
                            }
                            case kRemoveLifetime: {
                                auto rows = currentLifetimes();
                                const long k = std::lround(removeIndex_->GetNumber());
                                if (k >= 1 && k <= static_cast<long>(rows.size()) && rows.size() > 1) {
                                    rows.erase(rows.begin() + (k - 1));
                                    rebuildLifetimeRows(rows);
                                }
                                break;
                            }
                            case kFit: runFit(); break;
                            case kSaveResults: saveResults(); break;
                            case kUnzoom:
                                if (display_) display_->GetXaxis()->UnZoom();
                                redraw();
                                break;
                            default: break;
                        }
                        break;
                    case kCM_RADIOBUTTON:
                        if (parm1 >= kModeZoom && parm1 <= kModeBackground)
                            setMouseMode(static_cast<MouseMode>(parm1 - kModeZoom));
                        break;
                    case kCM_CHECKBUTTON:
                        if (parm1 == kLogY || parm1 == kShowComponents) redraw();
                        else if (parm1 == kUseFitStart || parm1 == kUseBgRange) applyAndEstimate();
                        break;
                    default: break;
                }
                break;
            case kC_TEXTENTRY:
                if (GET_SUBMSG(msg) == kTE_ENTER && parm1 >= kFitStart && parm1 <= kBgTo) applyAndEstimate();
                break;
            default: break;
        }
    } catch (const std::exception& e) {
        error(e.what());
    }
    return kTRUE;
}

} // namespace pals
