void ImageProcessingGray(unsigned char* pOut
	, unsigned char* pIn
	, size_t nWidth
	, size_t nHeight)
{
	for (int i = 0; i < nHeight; i++) {
		for (int j = 0; j < nWidth; j++) {
			if (pIn[i * nWidth + j] >= 123) {
				pOut[i * nWidth + j] = 0;
			}
			else {
				pOut[i * nWidth + j] = 255;
			}
		}
	}
}



void ImageProcessingGray(unsigned char* pOut
	, unsigned char* pIn
	, size_t nWidth
	, size_t nHeight)
{
	for (int i = 0; i < nHeight; i++) {
		for (int j = 0; j < nWidth; j++) {
			if (pIn[i * nWidth + j] >= 123) {
				pIn[i * nWidth + j] = 0;
			}
			else {
				pIn[i * nWidth + j] = 1;
			}
		}
	}

	int label = 1;
	int svaz = 4; // определение типа связи, ВАЖНО!!!
	std::map<int, int> remarks; // словарь, чтобы помещать в него значения для ремаркировки
	if (svaz == 4) {
		for (int i = 0; i < nHeight; i++) {
			for (int j = 0; j < nWidth; j++) {
				if (pIn[i * nWidth + j] != 0) {
					if (i > 0 && j > 0) {
						if (pIn[(i - 1) * nWidth + j] != 0 && pIn[i * nWidth + j - 1] != 0) {
							int tmp1 = std::min(pIn[(i - 1) * nWidth + j], pIn[i * nWidth + j - 1]);
							int tmp2 = std::max(pIn[(i - 1) * nWidth + j], pIn[i * nWidth + j - 1]);
							pIn[i * nWidth + j] = tmp1;
							if (tmp1 != tmp2) {
								remarks[tmp2] = tmp1; // в словарь по ключу добавляются элементы, которые надо заменить на ключ
							}

						}
						else if (pIn[(i - 1) * nWidth + j] != 0) {
							pIn[i * nWidth + j] = pIn[(i - 1) * nWidth + j];
						}
						else if (pIn[i * nWidth + j - 1] != 0) {
							pIn[i * nWidth + j] = pIn[i * nWidth + j - 1];
						}
						else {
							pIn[i * nWidth + j] = label++;
						}
					}
					else if (i > 0) {
						if (pIn[(i - 1) * nWidth + j] != 0) {
							pIn[i * nWidth + j] = pIn[(i - 1) * nWidth + j];
						}
						else {
							pIn[i * nWidth + j] = label++;
						}
					}
					else if (j > 0) {
						if (pIn[i * nWidth + j - 1] != 0) {
							pIn[i * nWidth + j] = pIn[i * nWidth + j - 1];
						}
						else {
							pIn[i * nWidth + j] = label++;
						}
					}
				}
			}
		}
		for (int i = label - 1; i > 0; i--) {
			if (remarks.count(i)) {
				while (remarks.count(remarks[i])) {
					int tmp = remarks[i];
					remarks[i] = remarks[tmp];
				}
			}
		}
		int levels = 255 / (label - remarks.size() - 1);	
		std::cout << "Number of objects: " << (label - remarks.size() - 1);
		for (int i = 0; i < nHeight; i++) {
			for (int j = 0; j < nWidth; j++) {
				if (remarks.count(pIn[i * nWidth + j])) {
					pOut[i * nWidth + j] = levels * remarks[pIn[i * nWidth + j]];
				}
				else {
					pOut[i * nWidth + j] = levels * pIn[i * nWidth + j];
				}
			}
		}
	}
	else if (svaz == 8) {
		for (int i = 0; i < nHeight; i++) {
			for (int j = 0; j < nWidth; j++) {
				if (pIn[i * nWidth + j] != 0) {
					if (i > 0 && j > 0) {
						int a = pIn[(i - 1) * nWidth + j];
						int b = pIn[i * nWidth + j - 1];
						int c = pIn[(i - 1) * nWidth + j - 1];

						if (a != 0 && b != 0 && c != 0) {
							int tmp1 = std::min(a, std::min(b, c));
							int tmp2 = std::max(a, std::max(b, c));
							int tmp3;
							if (tmp1 == a && tmp2 == b) {
								tmp3 = c;
							}
							else if (tmp1 == b && tmp2 == c) {
								tmp3 = a;
							}
							else {
								tmp3 = b;
							}
							pIn[i * nWidth + j] = tmp1;
							if (tmp1 != tmp2) {
								remarks[tmp2] = tmp1; // в словарь по ключу добавляются элементы, которые надо заменить на ключ
							}
							if (tmp1 != tmp3) {
								remarks[tmp3] = tmp1;
							}
						}
						else if (a != 0 && b != 0) {
							int tmp1 = std::min(a, b);
							int tmp2 = std::max(a, b);
							pIn[i * nWidth + j] = tmp1;
							if (tmp1 != tmp2) {
								remarks[tmp2] = tmp1;
							}
						}
						else if (b != 0 && c != 0) {
							int tmp1 = std::min(c, b);
							int tmp2 = std::max(c, b);
							pIn[i * nWidth + j] = tmp1;
							if (tmp1 != tmp2) {
								remarks[tmp2] = tmp1;
							}
						}
						else if (a != 0 && c != 0) {
							int tmp1 = std::min(c, a);
							int tmp2 = std::max(c, a);
							pIn[i * nWidth + j] = tmp1;
							if (tmp1 != tmp2) {
								remarks[tmp2] = tmp1;
							}
						}
						else if (a) {
							pIn[i * nWidth + j] = a;
						}
						else if (b) {
							pIn[i * nWidth + j] = b;
						}
						else if (c) {
							pIn[i * nWidth + j] = c;
						}
						else {
							pIn[i * nWidth + j] = label++;
						}
					}
					else if (i > 0) {
						if (pIn[(i - 1) * nWidth + j] != 0) {
							pIn[i * nWidth + j] = pIn[(i - 1) * nWidth + j];
						}
						else {
							pIn[i * nWidth + j] = label++;
						}
					}
					else if (j > 0) {
						if (pIn[i * nWidth + j - 1] != 0) {
							pIn[i * nWidth + j] = pIn[i * nWidth + j - 1];
						}
						else {
							pIn[i * nWidth + j] = label++;
						}
					}
				}
			}
		}
		for (int i = label - 1; i > 0; i--) {
			if (remarks.count(i)) {
				while (remarks.count(remarks[i])) {
					int tmp = remarks[i];
					remarks[i] = remarks[tmp];
				}
			}
		}
		int levels = 255 / (label - remarks.size() - 1);
		std::cout << "Number of objects: " << (label - remarks.size() - 1);
		for (size_t i = 0; i < nHeight; i++) {
			for (size_t j = 0; j < nWidth; j++) {
				if (remarks.count(pIn[i * nWidth + j])) {
					pOut[i * nWidth + j] = levels * remarks[pIn[i * nWidth + j]];
				}
				else {
					pOut[i * nWidth + j] = levels * pIn[i * nWidth + j];
				}
			}
		}
	}
	return;
}