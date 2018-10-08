//методы алгоритма
void simplex_method(vector<double> & matrix_x, std::string & string) {
	vector<vector<double>> matrix_a;
	vector<double> matrix_b;
	vector<double> matrix_c;
	bool limit, equality;

	parse_input(string, limit, equality, matrix_a, matrix_b, matrix_c);

	//make canon matrix
	if (!limit) {
		for (size_t i = 0; i < matrix_c.size(); ++i) {
			matrix_c.at(i) *= -1;
		}
	}

	if (!equality) {
		for (size_t i = 0; i < matrix_b.size(); ++i) {
			matrix_b.at(i) *= (-1);
			for (size_t j = 0; j < matrix_a.at(i).size(); ++j) {
				matrix_a.at(i).at(j) *= (-1);
			}
		}
	}

	vector<vector<double>> simplex_matrix = make_simplex_matrix(matrix_a, matrix_b, matrix_c);

	//two stages
	cout << "first stage:\n";
	if (!stages(simplex_matrix, 1)) {
		cout << "there are no feasible solutions\n";
	}

	print_simplex_matrix(simplex_matrix);
	cout << "\nsecond stage:\n";
	if (!stages(simplex_matrix, 2)) {
		cout << "there are no optimal solutions\n";
	}
	
	//inversion
	if (!limit) {
		for (size_t i = 0; i < simplex_matrix.at(simplex_matrix.size() - 1).size(); ++i) {
			simplex_matrix.at(simplex_matrix.size() - 1).at(i) *= -1;
		}
	}
	
	print_simplex_matrix(simplex_matrix);
	matrix_x = result(simplex_matrix);

	//finish check
	if (!check(matrix_a, matrix_b, matrix_c, matrix_x)) {
		throw std::invalid_argument("wrong result");
	}

	print_result(matrix_x);
}

bool stages(vector<vector<double>> & simplex_matrix, const size_t mode) {
	size_t row, column, supp_elem;
	while (!check_solution(simplex_matrix, mode, supp_elem)) {
		if (search_supporting_element(simplex_matrix, mode, supp_elem, row, column)) {
			do_Jordan_exceptions(simplex_matrix, row, column);
			print_simplex_matrix(simplex_matrix);
		}
		else {
			return false;
		}
	}
	return true;
}

bool check_solution(vector<vector<double>> & simplex_matrix, const size_t mode, size_t & supp_elem) {
	//check base solution
	if (mode == 1) {
		for (supp_elem = 1; supp_elem < simplex_matrix.size() - 1; ++supp_elem) {
			if (simplex_matrix.at(supp_elem).at(1) < 0) {
				return false;
			}
		}
	}
	//check optimal solution
	else if (mode == 2) {
		for (supp_elem = 1; supp_elem < simplex_matrix.at(simplex_matrix.size() - 1).size(); ++supp_elem) {
			if (simplex_matrix.at(simplex_matrix.size() - 1).at(supp_elem) > 0) {
				return false;
			}
		}
	}

	return true;
}

bool search_supporting_element(vector<vector<double>> & simplex_matrix, const size_t mode, size_t & supp_elem, size_t & row, size_t & column) {
	size_t counter = 0;
	if (mode == 1) {
		//поиск разрешающего стобца
		for (size_t k = 2; k < simplex_matrix.at(supp_elem).size(); ++k) {
			if (simplex_matrix.at(supp_elem).at(k) < 0) {
				column = k;
				break;
			}
			else {
				++counter;
			}
		}
		if (counter == simplex_matrix.at(supp_elem).size() - 2) {
			return false;
		}
		//поиск разрешающей строки
		double quotient = simplex_matrix.at(supp_elem).at(1) / simplex_matrix.at(supp_elem).at(column);
		row = supp_elem;

		for (size_t j = 1; j < simplex_matrix.size() - 1; ++j) {
			if (quotient > simplex_matrix.at(j).at(1) / simplex_matrix.at(j).at(column) &&
				simplex_matrix.at(j).at(1) / simplex_matrix.at(j).at(column) > 0) {
				quotient = simplex_matrix.at(j).at(1) / simplex_matrix.at(j).at(column);
				row = j;
			}
		}
	}
	else if (mode == 2) {
		//проверка столбца
		for (size_t k = 1; k < simplex_matrix.size() - 1; ++k) {
			if (simplex_matrix.at(k).at(supp_elem) > 0) {
				row = k;
				break;
			}
			else {
				++counter;
			}
		}
		if (counter == simplex_matrix.size() - 2) {
			return false;
		}

		//поиск разрешающей строки
		double quotient = simplex_matrix.at(row).at(1) / simplex_matrix.at(row).at(supp_elem);
		column = supp_elem;

		for (size_t j = 1; j < simplex_matrix.size() - 1; ++j) {
			if (quotient > simplex_matrix.at(j).at(1) / simplex_matrix.at(j).at(column) &&
				simplex_matrix.at(j).at(column) > 0) {
				quotient = simplex_matrix.at(j).at(1) / simplex_matrix.at(j).at(column);
				row = j;
			}
		}
	}
	return true;
}

void do_Jordan_exceptions(vector<vector<double>> & simplex_matrix, size_t row, size_t column) {
	//меняем имена переменных
	std::swap(simplex_matrix.at(row).at(0), simplex_matrix.at(0).at(column));
	//меняем таблицу
	for (size_t i = 1; i < simplex_matrix.size(); ++i) {
		if (i == row) {
			continue;
		}
		for (size_t j = 1; j < simplex_matrix.at(i).size(); ++j) {
			if (j == column) {
				continue;
			}
			simplex_matrix.at(i).at(j) -= (simplex_matrix.at(i).at(column) * simplex_matrix.at(row).at(j)) / simplex_matrix.at(row).at(column);

		}
	}

	for (size_t j = 1; j < simplex_matrix.at(row).size(); ++j) {
		if (j == column) {
			continue;
		}
		simplex_matrix.at(row).at(j) /= simplex_matrix.at(row).at(column);
	}

	for (size_t i = 1; i < simplex_matrix.size(); ++i) {
		if (i == row) {
			continue;
		}
		simplex_matrix.at(i).at(column) /= ((-1) * simplex_matrix.at(row).at(column));
	}

	simplex_matrix.at(row).at(column) = 1 / simplex_matrix.at(row).at(column);
}

bool check(vector<vector<double>> & matrix_a, vector<double> & matrix_b, vector<double> & matrix_c, vector<double> & matrix_x) {
	for (size_t i = 0; i < matrix_a.size(); ++i) {
		double sum = 0;
		for (size_t j = 0; j < matrix_x.size(); ++j) {
			sum += matrix_a.at(i).at(j) * matrix_x.at(j);
		}
		if (sum > matrix_b.at(i)) {
			return false;
		}
	//	cout << sum << '-' << matrix_b.at(i) << endl;
	}

	return true;
}

//методы создания
vector<vector<double>> make_simplex_matrix(vector<vector<double>> & matrix_a, vector<double> & matrix_b, vector<double> & matrix_c) {
	vector<vector<double>> simplex_matrix;

	for (size_t i = 0; i < matrix_b.size() + 2; ++i) {
		vector<double> mid_result;

		for (size_t j = 0; j < matrix_c.size() + 2; ++j) {
			if (i == 0 && j == 0) {
				mid_result.push_back(0);
			}
			else if (i == 0 && j > 0) {
				mid_result.push_back(j - 1);
			}
			else if (i > 0 && j == 0) {
				mid_result.push_back(matrix_c.size() + i);
			}
			else if (i > 0 && i < matrix_b.size() + 1 && j == 1) {
				mid_result.push_back(matrix_b.at(i - 1));
			}
			else if (i < matrix_a.size() + 1 && j > 1) {
				mid_result.push_back(matrix_a.at(i - 1).at(j - 2));
			}
			else if (i == matrix_a.size() + 1 && j == 1) {
				mid_result.push_back(0);
			}
			else if (i == matrix_a.size() + 1 && j > 1) {
				mid_result.push_back(matrix_c.at(j - 2) * (-1));
			}
		}
		simplex_matrix.push_back(mid_result);
	}

	return simplex_matrix;
}

vector<double> result(vector<vector<double>> & simplex_matrix) {
	size_t size = simplex_matrix.at(0).size() - 2;
	vector<double> matrix_result;

	for (size_t i = 0; i < size; ++i) {
		matrix_result.push_back(0);
	}

	for (size_t i = 1; i < simplex_matrix.size() - 1; ++i) {
		if (simplex_matrix.at(i).at(0) <= size) {
			int a = (int)simplex_matrix.at(i).at(0) - 1;
			matrix_result.at(a) = simplex_matrix.at(i).at(1);
		}
	}

	return matrix_result;
}

//методы вывода
void print_simplex_matrix(vector<vector<double>> & matrix) {
	std::ofstream fout("output.txt", std::ios::app);
	if (!fout.is_open()) {
		throw std::invalid_argument("wrong path");
	}
	cout << endl;
	for (size_t i = 0; i < matrix.size(); ++i) {
		for (size_t j = 0; j < matrix.at(i).size(); ++j) {
			if (i == 0 && j == 0) {
				cout << ' ' << '\t';
			}
				else if (i == 0 && j == 1) {
				cout << " S" << '\t';
			}
			else if (i == 0 && j > 1 || j == 0 && i > 0 && i < matrix.size() - 1) {
				cout << "x" << matrix.at(i).at(j) << '\t';
			}
			else if (i == matrix.size() - 1 && j == 0) {
				cout << "F" << '\t';
			}
			else {
				if (matrix.at(i).at(j) >= 0) {
					cout << ' ';
					fout << ' ';
				}
				if (matrix.at(i).at(j) == 0) {
					cout << 0 << '\t';
					fout << 0 << '\t';
				}
				else {
					cout << matrix.at(i).at(j) << '\t';
					fout << matrix.at(i).at(j) << '\t';
				}
			}
			if (i == 0 && j == matrix.at(i).size() - 1) {
				cout << endl;
				fout << endl;
			}
		}
		cout << endl;
		fout << endl;
	}

	cout << "F = " << matrix.at(matrix.size() - 1).at(1) << endl;
	fout << "F = " << matrix.at(matrix.size() - 1).at(1) << endl;

	fout.close();
}

void print_result(vector<double> & matrix_x) {
	std::ofstream fout("output.txt", std::ios::app);
	if (!fout.is_open()) {
		throw std::invalid_argument("wrong path");
	}

	for (auto i : matrix_x) {
		cout << i << '\t';
		fout << i << '\t';
	}
	cout << endl;
	fout << endl;

	fout.close();
}

//методы считывания
void parse_input(std::string & string, bool & limit, bool & equality, vector<vector<double>> & matrix_a, vector<double> & matrix_b, vector<double> & matrix_c) {
	std::istringstream stream(string);

	limit = parse_limit(stream);
	equality = parse_equality(stream);

	vector<bool> counter = { true, true, true };

	for (size_t i = 0; i < 3; ++i) {
		std::string new_string;
		getline(stream, new_string, '=');

		if (new_string == "c") {
			if (counter.at(0) == false) {
				throw std::logic_error("wrong string");
			}
			matrix_c = parse_array(stream);
			counter.at(0) = false;
		}
		else if (new_string == "A") {
			if (counter.at(1) == false) {
				throw std::logic_error("wrong string");
			}
			matrix_a = parse_double_array(stream);
			counter.at(1) = false;
		}
		else if (new_string == "b") {
			if (counter.at(2) == false) {
				throw std::logic_error("wrong string");
			}
			matrix_b = parse_array(stream);
			counter.at(2) = false;
		}
		else {
			throw std::logic_error("wrong string");
		}
	}
}

bool parse_limit(std::istringstream & stream) {
	std::string new_string;
	std::getline(stream, new_string, ';');

	if (new_string == "min") {
		return true;
	}
	else if (new_string == "max") {
		return false;
	}
	else {
		throw std::logic_error("wrong string");
	}

	return false;
}

bool parse_equality(std::istringstream & stream) {
	std::string new_string;
	std::getline(stream, new_string, ';');

	if (new_string == "<=") {
		return true;
	}
	else if (new_string == ">=") {
		return false;
	}
	else {
		throw std::logic_error("wrong string");
	}

	return false;
}

vector<double> parse_array(std::istringstream & stream) {
	vector<double> matrix;
	char symbol;
	if (stream >> symbol && symbol == '{') {
		double perem;

		while (stream >> perem) {
			matrix.push_back(perem);

			if (stream >> symbol && symbol == ',') {
				continue;
			}
			else if (symbol == '}') {
				break;
			}
			else {
				throw std::logic_error("wrong string");
			}
		}
	}
	else {
		throw std::logic_error("wrong string");
	}

	return matrix;
}

vector<vector<double>> parse_double_array(std::istringstream & stream) {
	vector<vector<double>> matrix;
	char symbol;

	if (stream >> symbol && symbol == '{') {
		while (symbol != '}') {
			matrix.push_back(parse_array(stream));

			if (stream >> symbol && symbol == ',') {
				continue;
			}
			else if (symbol == '}') {
				break;
			}
			else {
				throw std::logic_error("wrong string");
			}
		}
	}
	else {
		throw std::logic_error("wrong string");
	}

	return matrix;
}

std::string parse_file(const std::string & path_to_file) {
	std::ifstream stream(path_to_file);

	if (!stream.is_open()) {
		throw std::invalid_argument("wrong path");
	}

	std::string string;
	char symbol;
	while (stream >> symbol) {
		string += symbol;
	}

	stream.close();

	return string;
}
