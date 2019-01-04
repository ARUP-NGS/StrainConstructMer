// Type definitions for ag-grid-community v19.1.1
// Project: http://www.ag-grid.com/
// Definitions by: Niall Crosby <https://github.com/ag-grid/>
import { Column } from "../entities/column";
import { RowNode } from "../entities/rowNode";
import { GridCell } from "../entities/gridCell";
import { ICellEditorComp } from "./cellEditors/iCellEditor";
import { Component } from "../widgets/component";
import { ICellRendererComp } from "./cellRenderers/iCellRenderer";
import { Beans } from "./beans";
import { RowComp } from "./rowComp";
export declare class CellComp extends Component {
    static DOM_DATA_KEY_CELL_COMP: string;
    private eCellWrapper;
    private eParentOfValue;
    private beans;
    private column;
    private rowNode;
    private eParentRow;
    private gridCell;
    private rangeCount;
    private usingWrapper;
    private includeSelectionComponent;
    private includeRowDraggingComponent;
    private cellFocused;
    private editingCell;
    private cellEditorInPopup;
    private hideEditorPopup;
    private lastIPadMouseClickEvent;
    private usingCellRenderer;
    private cellRendererType;
    private cellRenderer;
    private cellRendererGui;
    private cellEditor;
    private autoHeightCell;
    private firstRightPinned;
    private lastLeftPinned;
    private rowComp;
    private rangeSelectionEnabled;
    private value;
    private valueFormatted;
    private colsSpanning;
    private rowSpan;
    private tooltip;
    private scope;
    private readonly printLayout;
    private cellEditorVersion;
    private cellRendererVersion;
    constructor(scope: any, beans: Beans, column: Column, rowNode: RowNode, rowComp: RowComp, autoHeightCell: boolean, printLayout: boolean);
    getCreateTemplate(): string;
    private getStylesForRowSpanning;
    afterAttached(): void;
    private onColumnHover;
    private onCellChanged;
    private getCellLeft;
    private getCellWidth;
    private onFlashCells;
    private setupColSpan;
    private getColSpanningList;
    private onDisplayColumnsChanged;
    private getInitialCssClasses;
    getInitialValueToRender(): string;
    getRenderedRow(): RowComp;
    isSuppressNavigable(): boolean;
    getCellRenderer(): ICellRendererComp;
    getCellEditor(): ICellEditorComp;
    refreshCell(params?: {
        suppressFlash?: boolean;
        newData?: boolean;
        forceRefresh?: boolean;
    }): void;
    flashCell(): void;
    private animateCell;
    private replaceContentsAfterRefresh;
    private angular1Compile;
    private postProcessStylesFromColDef;
    private preProcessStylesFromColDef;
    private processStylesFromColDef;
    private postProcessClassesFromColDef;
    private preProcessClassesFromColDef;
    private processClassesFromColDef;
    private putDataIntoCellAfterRefresh;
    attemptCellRendererRefresh(): boolean;
    private refreshToolTip;
    private valuesAreEqual;
    private getToolTip;
    private processCellClassRules;
    private postProcessCellClassRules;
    private preProcessCellClassRules;
    setUsingWrapper(): void;
    private chooseCellRenderer;
    private createCellRendererInstance;
    private afterCellRendererCreated;
    private attachCellRenderer;
    private createCellRendererParams;
    private formatValue;
    private getValueToUse;
    private getValueAndFormat;
    private getValue;
    onMouseEvent(eventName: string, mouseEvent: MouseEvent): void;
    dispatchCellContextMenuEvent(event: Event): void;
    private createEvent;
    private onMouseOut;
    private onMouseOver;
    private onCellDoubleClicked;
    startRowOrCellEdit(keyPress?: number, charPress?: string): void;
    isCellEditable(): boolean;
    startEditingIfEnabled(keyPress?: number, charPress?: string, cellStartedEdit?: boolean): void;
    private afterCellEditorCreated;
    private addInCellEditor;
    private addPopupCellEditor;
    private onPopupEditorClosed;
    private setInlineEditingClass;
    private createCellEditorParams;
    private stopEditingAndFocus;
    private parseValue;
    focusCell(forceBrowserFocus?: boolean): void;
    setFocusInOnEditor(): void;
    isEditing(): boolean;
    onKeyDown(event: KeyboardEvent): void;
    doesUserWantToCancelKeyboardEvent(event: KeyboardEvent): boolean;
    setFocusOutOnEditor(): void;
    private onNavigationKeyPressed;
    private onShiftRangeSelect;
    private onTabKeyDown;
    private onBackspaceOrDeleteKeyPressed;
    private onEnterKeyDown;
    private navigateAfterEdit;
    private onF2KeyDown;
    private onEscapeKeyDown;
    onKeyPress(event: KeyboardEvent): void;
    private onSpaceKeyPressed;
    private onMouseDown;
    private isDoubleClickOnIPad;
    private onCellClicked;
    private doIeFocusHack;
    private createGridCellVo;
    getGridCell(): GridCell;
    getParentRow(): HTMLElement;
    setParentRow(eParentRow: HTMLElement): void;
    getColumn(): Column;
    detach(): void;
    destroy(): void;
    private onLeftChanged;
    private modifyLeftForPrintLayout;
    private onWidthChanged;
    private getRangeClasses;
    private onRowIndexChanged;
    private onRangeSelectionChanged;
    private onFirstRightPinnedChanged;
    private onLastLeftPinnedChanged;
    private populateTemplate;
    private addRowDragging;
    private addSelectionCheckbox;
    private addDomData;
    private onCellFocused;
    stopRowOrCellEdit(cancel?: boolean): void;
    stopEditing(cancel?: boolean): void;
}
//# sourceMappingURL=cellComp.d.ts.map