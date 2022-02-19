#include <windows.h>
#include <d2d1.h>
#pragma comment(lib, "d2d1")

#include "FluidSim.h"
#include "Fluid.h"

#define SHOW_VEL false


template <class T> void SafeRelease(T** ppT)
{
    if (*ppT)
    {
        (*ppT)->Release();
        *ppT = NULL;
    }
}

class MainWindow : public BaseWindow<MainWindow>
{
    ID2D1Factory* pFactory;
    ID2D1HwndRenderTarget* pRenderTarget;
    ID2D1SolidColorBrush* pBrush;
    //D2D1_ELLIPSE            ellipse;

    void    CalculateLayout();
    HRESULT CreateGraphicsResources();
    void    DiscardGraphicsResources();
    void    OnPaint();
    void    Resize();

    POINT _prevPos;

    // other Variables
    Fluid _fluid;

    float _theta;

public:

    MainWindow() : pFactory(NULL), pRenderTarget(NULL), pBrush(NULL)
    {
        _fluid = Fluid(winSize / gridSize, gridSize, 0, 0, 5000, 0, 0.1);

        GetCursorPos(&_prevPos);

        _theta = 0;

        for (int i = 0; i < winSize / gridSize; i++)
        {
            for (int j = 0; j < winSize / gridSize; j++)
            {
                if (i > 10 && i < 15 && j > 10 && j < 15)
                {
                    _fluid.addPrevDensity(i, j, 10);
                }

                if (i > 5 && i < 20 && j > 2 && j < 20)
                {
                    _fluid.setPreviousVelocity(i, j, 10, 10);
                }
            }
        }
    }

    PCWSTR  ClassName() const { return L"Circle Window Class"; }
    LRESULT HandleMessage(UINT uMsg, WPARAM wParam, LPARAM lParam);
};

// Recalculate drawing layout when the size of the window changes.

/*void MainWindow::CalculateLayout()
{
    if (pRenderTarget != NULL)
    {
        D2D1_SIZE_F size = pRenderTarget->GetSize();
        const float x = size.width / 2;
        const float y = size.height / 2;
        const float radius = min(x, y);
        ellipse = D2D1::Ellipse(D2D1::Point2F(x, y), radius, radius);
    }
}*/

HRESULT MainWindow::CreateGraphicsResources()
{
    HRESULT hr = S_OK;
    if (pRenderTarget == NULL)
    {
        RECT rc;
        GetClientRect(m_hwnd, &rc);

        D2D1_SIZE_U s = D2D1::SizeU(rc.right, rc.bottom);

        hr = pFactory->CreateHwndRenderTarget(
            D2D1::RenderTargetProperties(),
            D2D1::HwndRenderTargetProperties(m_hwnd, s),
            &pRenderTarget);

        if (SUCCEEDED(hr))
        {
            const D2D1_COLOR_F color = D2D1::ColorF(1.0f, 1.0f, 0);
            hr = pRenderTarget->CreateSolidColorBrush(color, &pBrush);

            if (SUCCEEDED(hr))
            {
            }
        }
    }
    return hr;
}

void MainWindow::DiscardGraphicsResources()
{
    SafeRelease(&pRenderTarget);
    SafeRelease(&pBrush);
}

void MainWindow::OnPaint()
{
    HRESULT hr = CreateGraphicsResources();
    if (SUCCEEDED(hr))
    {
        PAINTSTRUCT ps;
        BeginPaint(m_hwnd, &ps);

        pRenderTarget->BeginDraw();

        //pRenderTarget->Clear(D2D1::ColorF(D2D1::ColorF::Black));

        //paint speed
        for (int i = 0; i < winSize / gridSize; i++)
        {
            for (int j = 0; j < winSize / gridSize; j++)
            {
                D2D1_RECT_F r = { gridSize * i, gridSize * j, gridSize * (i + 5), gridSize * (j + 5) };

                float _d = _fluid.getDensity(i, j);
                float _vel = sqrt(_fluid.getVelocity(i, j, 0) * _fluid.getVelocity(i, j, 0) + _fluid.getVelocity(i, j, 1) * _fluid.getVelocity(i, j, 1));
                
                _fluid.addDensity(i, j, _fluid.getDensity(i, j) - 0.001);

                pBrush->SetColor(D2D1::ColorF(_d, _d, 0));

                pRenderTarget->FillRectangle(r, pBrush);

                if (SHOW_VEL)
                {
                    pBrush->SetColor(D2D1::ColorF(D2D1::ColorF::White));

                    D2D1_POINT_2F _start = { i * gridSize, j * gridSize };
                    D2D1_POINT_2F _end = { (i + _fluid.getVelocity(i, j, 0) * 0.1) * gridSize , (j + _fluid.getVelocity(i, j, 1) * 0.1) * gridSize };

                    pRenderTarget->DrawLine(_start, _end, pBrush);
                    //_fluid.setVelocity(i, j, _fluid.getVelocity(i, j, 0) + 0.01, _fluid.getVelocity(i, j, 1) + 0.01);
                }
            }
        }
        POINT _pos;

        GetCursorPos(&_pos);

        RECT r;

        GetWindowRect(m_hwnd, &r);

        if (_pos.x > r.left + 2 * gridSize
            && _pos.x < r.right - 2 * gridSize
            && _pos.y > r.top + 50
            && _pos.y < r.bottom - 2 * gridSize && (GetKeyState(VK_LBUTTON) & 0x100) != 0)
        {
            int _x = (_pos.x - r.left - 2 * gridSize) / gridSize;
            int _y = (_pos.y - r.top - 50) / gridSize;

           for (int i = _x; i < _x + 4; i++)
                for (int j = _y; j < _y + 4; j++)
                {
                    _fluid.addDensity(i, j, 3);
                    _fluid.setVelocity(i, j, 50, 50);
                }
        }

        _fluid.simulate(0.1);

        pBrush->SetColor(D2D1::ColorF(D2D1::ColorF::White));

        int c = winSize / (2 * gridSize);
        int rad = 20;

        /*for (int x = c - rad; x < c + rad; x++)
            for (int y = c - rad; y < c + rad; y++)
            {
                float _norm = sqrt((float) (x - c) * (float)(x - c) + (float) (y - c) * (float)(y - c));

                if (_norm <= 3 && _norm >= 1)
                {
                    D2D1_RECT_F rect = { gridSize * x, gridSize * y, gridSize * (x + 5), gridSize * (y + 5) };
                    pRenderTarget->FillRectangle(rect, pBrush);
                }
            }*/

        hr = pRenderTarget->EndDraw();

        if (FAILED(hr) || hr == D2DERR_RECREATE_TARGET)
        {
            DiscardGraphicsResources();
        }

        EndPaint(m_hwnd, &ps);
    }
}

void MainWindow::Resize()
{
    if (pRenderTarget != NULL)
    {
        RECT rc;
        GetClientRect(m_hwnd, &rc);

        D2D1_SIZE_U s = D2D1::SizeU(rc.right, rc.bottom);

        pRenderTarget->Resize(s);
       // CalculateLayout();
        InvalidateRect(m_hwnd, NULL, FALSE);
    }
}

int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE, PWSTR, int nCmdShow)
{
    MainWindow win;

    if (!win.Create(L"Circle", WS_OVERLAPPEDWINDOW))
    {
        return 0;
    }

    ShowWindow(win.Window(), nCmdShow);

    // Run the message loop.

    MSG msg = { };
    while (GetMessage(&msg, NULL, 0, 0))
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    return 0;
}

LRESULT MainWindow::HandleMessage(UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    switch (uMsg)
    {
    case WM_CREATE:
        SetTimer(m_hwnd, 0, 50, 0);

        if (FAILED(D2D1CreateFactory(
            D2D1_FACTORY_TYPE_SINGLE_THREADED, &pFactory)))
        {
            return -1;  // Fail CreateWindowEx.
        }
        return 0;

    case WM_DESTROY:
        DiscardGraphicsResources();
        SafeRelease(&pFactory);
        PostQuitMessage(0);
        return 0;

    case WM_PAINT:
        OnPaint();
        return 0;

    case WM_TIMER:
        OnPaint();
        return 0;

    case WM_SIZE:
        Resize();
        return 0;
    }
    return DefWindowProc(m_hwnd, uMsg, wParam, lParam);
}