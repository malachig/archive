#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <stdlib.h>
#include <stdio.h>

#define IE_EXPLORER "IEFrame"
#define IE_STATUSBAR "msctls_statusbar32"
#define MAX_URL_BUFFER 2084
#define DELAY 50

enum {
	HWND_IE,
		HWND_MSCTLS,
		HWND_SIZE_OF
	};

void error(char* s, DWORD dwCode)
{
	printf("%s", s);
	exit(dwCode);
}

void info()
{
	printf("IE_Spy by ___LiquidBinary___\n");
	printf("liquidbinary@linuxmail.org\n");
	printf("<CTRL-C> to quit\n\n");
}

int main(void)
{
	HWND hWnds[HWND_SIZE_OF];
	char sBuffer[MAX_URL_BUFFER];
	char sURL[MAX_URL_BUFFER];

	info();
	hWnds[HWND_IE]=FindWindow(IE_EXPLORER, NULL);
	if (!hWnds[HWND_IE])
		error("IE not opened ...\n", 0);

	hWnds[HWND_MSCTLS]=FindWindowEx(hWnds[HWND_IE], NULL, IE_STATUSBAR, NULL);

	if(!hWnds[HWND_MSCTLS])
		error("Cannot locate status bar ...\n", 0);

	printf("Logging all IE URL requests ...\n");
	SendMessage(hWnds[HWND_MSCTLS], WM_GETTEXT, MAX_URL_BUFFER, (LPARAM) sBuffer);

	printf("%s\n", sBuffer);
		for( ;; )
		{
			SendMessage(hWnds[HWND_MSCTLS], WM_GETTEXT, MAX_URL_BUFFER, (LPARAM) sURL);
			Sleep(DELAY);
			if(lstrcmp(sURL,sBuffer)==0)
				continue;
			else
			{
				printf("%s\n", sURL);
				lstrcpy(sBuffer, sURL);
			}
		}
		return 0;
}
