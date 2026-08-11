function configureGraphics()
    %% Latex
    set(groot, 'defaultAxesTickLabelInterpreter','latex');  
    set(groot, 'defaultBubblelegendInterpreter','latex');  
    set(groot, 'defaultColorbarTickLabelInterpreter','latex');
    set(groot, 'defaultConstantlineInterpreter','latex');  
    set(groot, 'defaultGraphplotInterpreter','latex');  
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultPolaraxesTickLabelInterpreter','latex');  
    set(groot, 'defaultTextInterpreter','latex');  
    set(groot, 'defaultTextarrowshapeInterpreter','latex');
    set(groot, 'defaultTextboxshapeInterpreter','latex');

    %% Font size
    DEFAULT_FONT_SIZE = 18;
    set(groot,'defaultAxesFontSize',DEFAULT_FONT_SIZE)
    set(groot,'defaultBubblelegendFontSize',DEFAULT_FONT_SIZE)
    set(groot,'defaultColorbarFontSize',DEFAULT_FONT_SIZE)
    set(groot,'defaultConstantlineFontSize',DEFAULT_FONT_SIZE)
    set(groot,'defaultGeoaxesFontSize',DEFAULT_FONT_SIZE)
    set(groot,'defaultGraphplotEdgeFontSize',DEFAULT_FONT_SIZE)
    set(groot,'defaultGraphplotNodeFontSize',DEFAULT_FONT_SIZE)
    set(groot,'defaultLegendFontSize',DEFAULT_FONT_SIZE)
    set(groot,'defaultPolaraxesFontSize',DEFAULT_FONT_SIZE)
    set(groot,'defaultTextFontSize',DEFAULT_FONT_SIZE)
    set(groot,'defaultTextarrowshapeFontSize',DEFAULT_FONT_SIZE)
    set(groot,'defaultTextboxshapeFontSize',DEFAULT_FONT_SIZE)
    set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.0);
end